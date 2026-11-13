# **************************************************************************
# *
# * Authors:    Verónica Gamo
# *             Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *             Blanca Pueche (blanca.pueche@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import os
import json
import re
import webbrowser
import subprocess

import matplotlib.pyplot as plt

import pyworkflow.viewer as pwviewer
import pyworkflow.utils as pwutils
from pyworkflow.protocol import params

from pwchem import Plugin as pwchemPlugin
from pwchem.constants import OPENBABEL_DIC

from retrosynthesis.protocols.protocol_retrosynthesis import ProtChemAiZynthFinder


class ProtChemAiZynthFinderViewer(pwviewer.ProtocolViewer):
    """Viewer for AiZynthFinder retrosynthesis results."""

    _label = 'Retrosynthesis Viewer'
    _targets = [ProtChemAiZynthFinder]
    _environments = [pwviewer.DESKTOP_TKINTER]

    def _defineParams(self, form):

        form.addSection(label='Retrosynthesis routes')

        form.addParam(
            'selectedRoute',
            params.EnumParam,
            choices=self.getRouteChoices(),
            default=0,
            label='Select route:',
            help='Select a retrosynthesis route to inspect.'
        )

        form.addParam(
            'showRoute',
            params.LabelParam,
            label='Show reaction scheme:',
            help='Open the reaction scheme of the selected route.'
        )

        form.addParam(
            'showRouteInfo',
            params.LabelParam,
            label='Show route information:',
            help='Display statistics for the selected route.'
        )

        form.addParam(
            'showPrecursors',
            params.LabelParam,
            label='Open precursor molecules in PyMOL:',
            help='Visualize the terminal precursor molecules of the selected route.'
        )

        form.addParam(
            'showScores',
            params.LabelParam,
            label='Show route score distribution:',
            help='Display a histogram of the scores of all retrosynthesis routes.'
        )

    # ------------------------------------------------------------------
    # ROUTE INFORMATION
    # ------------------------------------------------------------------

    def getRoutesSummary(self):

        summary_file = self.protocol._getExtraPath(
            'routes_summary.json'
        )

        if not os.path.exists(summary_file):
            return []

        with open(summary_file, 'r') as f:
            return json.load(f)

    def getRouteChoices(self):

        routes = self.getRoutesSummary()

        if not routes:
            return ['No routes available']

        choices = []

        for route in routes:

            route_number = route['route']

            status = (
                'Solved' if route.get('solved', False)
                else 'Unsolved'
            )

            score = route.get('score', 'N/A')

            choices.append(
                'Route {} ({}, score: {})'.format(
                    route_number,
                    status,
                    score
                )
            )

        return choices

    def getSelectedRoute(self):

        routes = self.getRoutesSummary()

        if not routes:
            raise ValueError('No retrosynthesis routes are available.')

        selected_index = self.selectedRoute.get()

        if selected_index >= len(routes):
            raise ValueError('Invalid route selection.')

        return routes[selected_index]['route']

    def getSelectedRouteInfo(self):

        route_number = self.getSelectedRoute()

        for route in self.getRoutesSummary():

            if route['route'] == route_number:
                return route

        raise ValueError(
            'Information for Route {} was not found.'.format(
                route_number
            )
        )

    # ------------------------------------------------------------------
    # VISUALIZATION ACTIONS
    # ------------------------------------------------------------------

    def _getVisualizeDict(self):

        return {
            'showRoute': self._showRoute,
            'showRouteInfo': self._showRouteInfo,
            'showPrecursors': self._showPrecursors,
            'showScores': self._showScores
        }

    # ------------------------------------------------------------------
    # REACTION SCHEME
    # ------------------------------------------------------------------

    def _showRoute(self, paramName=None):

        route_number = self.getSelectedRoute()

        # Scipion Route 1 corresponds to route000.png.
        image_number = route_number - 1

        image_file = self.protocol._getExtraPath(
            'route{:03d}.png'.format(image_number)
        )

        if not os.path.exists(image_file):
            raise FileNotFoundError(
                'Reaction image not found: {}'.format(image_file)
            )

        subprocess.Popen(
            ['xdg-open', image_file]
        )
    # ------------------------------------------------------------------
    # ROUTE STATISTICS
    # ------------------------------------------------------------------

    def _showRouteInfo(self, paramName=None):
        import tkinter as tk
        from tkinter import ttk

        route = self.getSelectedRouteInfo()

        route_number = route['route']
        status = 'Solved' if route.get('solved', False) else 'Unsolved'

        info = [
            'Status: {}'.format(status),
            'Score: {}'.format(route.get('score', 'N/A')),
            'Reaction steps: {}'.format(route.get('reactions', 'N/A')),
            'Terminal precursors: {}'.format(route.get('precursors', 0)),
            'Precursors in stock: {}'.format(route.get('precursors_in_stock', 0)),
            'Precursors not in stock: {}'.format(route.get('precursors_not_in_stock', 0)),
            '',
            'Precursor SMILES:',
            ''
        ]

        info.extend(route.get('precursor_smiles', []))

        # Create popup window
        popup = tk.Toplevel()
        popup.title('Route {} - Information'.format(route_number))
        popup.geometry('750x500')
        popup.minsize(500, 300)

        # Main frame
        frame = ttk.Frame(popup, padding=10)
        frame.pack(fill='both', expand=True)

        # Route title
        title = ttk.Label(
            frame,
            text='Retrosynthesis Route {}'.format(route_number),
            font=('TkDefaultFont', 14, 'bold')
        )
        title.pack(anchor='w', pady=(0, 10))

        # Scrollable text area
        text_frame = ttk.Frame(frame)
        text_frame.pack(fill='both', expand=True)

        scrollbar = ttk.Scrollbar(text_frame)
        scrollbar.pack(side='right', fill='y')

        text_widget = tk.Text(
            text_frame,
            wrap='word',
            font=('TkFixedFont', 10),
            yscrollcommand=scrollbar.set,
            padx=8,
            pady=8
        )
        text_widget.pack(side='left', fill='both', expand=True)

        scrollbar.config(command=text_widget.yview)

        text_widget.insert('1.0', '\n'.join(info))
        text_widget.config(state='disabled')

        # Close button
        ttk.Button(
            frame,
            text='Close',
            command=popup.destroy
        ).pack(anchor='e', pady=(10, 0))

        popup.transient()
        popup.grab_set()

    # ------------------------------------------------------------------
    # PRECURSOR MOLECULES IN PYMOL
    # ------------------------------------------------------------------

    def _showPrecursors(self, paramName=None):

        route_number = self.getSelectedRoute()

        output_name = 'outputRoute{}'.format(route_number)

        if not hasattr(self.protocol, output_name):
            raise ValueError(
                'Output {} was not found.'.format(output_name)
            )

        output_set = getattr(self.protocol, output_name)

        if output_set is None or output_set.getSize() == 0:
            raise ValueError(
                'No precursor molecules were found for Route {}.'.format(
                    route_number
                )
            )

        # Create a directory for PyMOL scripts.
        pml_dir = self.protocol._getExtraPath('pmls')

        if not os.path.exists(pml_dir):
            os.makedirs(pml_dir)

        pml_file = os.path.join(
            pml_dir,
            'route{}_precursors.pml'.format(route_number)
        )

        pml_lines = []

        for idx, mol in enumerate(output_set, start=1):

            mol_file = os.path.abspath(mol.getFileName())

            if not os.path.exists(mol_file):
                print(
                    'Warning: precursor file not found: {}'.format(
                        mol_file
                    )
                )
                continue

            # Use a PyMOL-compatible object name.
            # Keep the original chemical name in Scipion.
            mol_name = mol.getMolName()

            safe_name = re.sub(
                r'[^a-zA-Z0-9_]', '_', str(mol_name)
            )

            safe_name = re.sub(r'_+', '_', safe_name).strip('_')

            if not safe_name:
                safe_name = 'Precursor{}'.format(idx)

            # Add the index to avoid duplicated object names.
            object_name = '{}_{}'.format(safe_name, idx)

            pml_lines.extend([
                'load "{}", {}'.format(mol_file, object_name),
                'hide everything, {}'.format(object_name),
                'show sticks, {}'.format(object_name)
            ])

        if not pml_lines:
            raise ValueError(
                'No valid precursor SDF files were found.'
            )

        pml_lines.append('zoom')

        with open(pml_file, 'w') as f:
            f.write('\n'.join(pml_lines))

        self.launchPyMOL(pml_file)

    def launchPyMOL(self, pml_file):

        pymol_bin = pwchemPlugin.getEnvPath(
            OPENBABEL_DIC,
            'bin/pymol'
        )

        if not os.path.exists(pymol_bin):
            raise FileNotFoundError(
                'PyMOL executable not found: {}'.format(pymol_bin)
            )

        environ = pwutils.Environ(os.environ)

        environ.set(
            'PATH',
            pwchemPlugin.getProgramHome(
                OPENBABEL_DIC,
                path='bin'
            ),
            position=pwutils.Environ.BEGIN
        )

        subprocess.Popen(
            [pymol_bin, os.path.abspath(pml_file)],
            cwd=os.path.dirname(pml_file),
            env=environ
        )

    # ------------------------------------------------------------------
    # ROUTE SCORE DISTRIBUTION
    # ------------------------------------------------------------------

    def _showScores(self, paramName=None):

        routes = self.getRoutesSummary()

        if not routes:
            raise ValueError('No retrosynthesis routes are available.')

        route_numbers = []
        scores = []

        for route in routes:

            score = route.get('score')

            if score is None:
                continue

            try:
                route_numbers.append(route['route'])
                scores.append(float(score))

            except (ValueError, TypeError):
                continue

        if not scores:
            raise ValueError('No valid route scores were found.')

        selected_route = self.getSelectedRoute()

        fig, ax = plt.subplots(figsize=(12, 6))
        bars = ax.bar(route_numbers, scores, edgecolor='black', alpha=0.8)

        for bar, route_number in zip(bars, route_numbers):
            if route_number == selected_route:
                bar.set_edgecolor('red')
                bar.set_linewidth(2.5)

        # Zoom the Y-axis around the observed score range
        min_score = min(scores)
        max_score = max(scores)
        score_range = max_score - min_score

        if score_range > 0:
            margin = score_range * 0.15
        else:
            margin = abs(max_score) * 0.05 or 0.01

        ax.set_ylim(min_score - margin, max_score + margin)

        ax.set_xlabel('Route number', fontsize=12)
        ax.set_ylabel('AiZynthFinder state score', fontsize=12)
        ax.set_title('Retrosynthesis route scores', fontsize=14)
        ax.set_xticks(route_numbers)
        ax.grid(axis='y', linestyle='--', alpha=0.4)
        ax.set_axisbelow(True)

        plt.tight_layout()
        plt.show()