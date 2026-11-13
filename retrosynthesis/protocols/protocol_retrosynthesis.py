################ VERÓNICA GAMO PAREJO ##############################

# General imports 
import gzip, json, os, shutil, csv, re
from urllib.request import urlopen

# Specific imports
from pyworkflow.protocol.params import PointerParam, STEPS_PARALLEL, StringParam
from pyworkflow.utils.path import moveFile, cleanPath
from pwem.protocols import EMProtocol

from pwchem.utils import *
from pwchem.constants import RDKIT_DIC
from pwchem.objects import SmallMolecule, SetOfSmallMolecules

from retrosynthesis.constants import AIZYNTHFINDER_DIC
from retrosynthesis import Plugin

class ProtChemAiZynthFinder(EMProtocol):

    """ Retrosynthesis of molecules with AiZynthFinder"""
    
    _label = 'aiZynthFinder'
    
    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL
        
    def _defineParams(self, form):
    
        form.addSection(label='Input')
                      
        form.addParam('inputSet', PointerParam, pointerClass="SetOfSmallMolecules",
                      label='Set of molecules to use:', allowsNull=False,
                      help='Select the set of small molecules for retrosynthesis.')
        
        form.addParam('inputLigand', StringParam, label='Selected just one ligand file:',
                      help='Relative path to the ligand file for the retrosynthesis.')
        
                   
    def _insertAllSteps(self):
        self._insertFunctionStep('identifyStep')
        self._insertFunctionStep('createOutputStep')

    def identifyStep(self):
        myLigand = str(self.inputLigand)
        ligand_file = None

        for mol in self.inputSet.get():
            if str(mol.getMolName()) in myLigand:
                ligand_file = mol.getFileName()
                break

        if not ligand_file:
            raise ValueError("No ligand file provided.")

        ligand_smiles = self.getSMI(ligand_file, 1)
        print("SMILES for retrosynthesis:", ligand_smiles)
        self.runAiZynthFinder(ligand_smiles)

    def getSMI(self, mol, nt):
        """Generate a SMILES representation from a molecule file."""
        fnSmall = os.path.abspath(mol)
        fnRoot, ext = os.path.splitext(os.path.basename(fnSmall))
        if ext == '.smi':
            fnOut = fnSmall
        else:
            outDir = os.path.abspath(self._getExtraPath())
            fnOut = os.path.abspath(self._getExtraPath(fnRoot + '.smi'))
            args = ' -i "{}" -of smi -o "{}" --outputDir "{}" -nt {}'.format(
                fnSmall, fnOut, outDir, nt
            )
            Plugin.runScript(
                self,
                'rdkit_IO.py',
                args,
                env=RDKIT_DIC,
                cwd=outDir
            )
        return self.parseSMI(fnOut)
        
    def parseSMI(self, smiFile):
        smi = None
        with open(smiFile) as f:
            for line in f:
                smi = line.split()[0].strip()
                if not smi.lower() == 'smiles':
                    break
        return smi
    
    def runAiZynthFinder(self, smi):
        # SAVE SMILE TO TXT FILE
        smiles_file = self._getExtraPath("smiles.txt")
        with open(smiles_file, 'w') as f:
            f.write(smi)
            print(f"SMILES {smi} saved to {smiles_file}")

        # RUN AIZYNTHFINDER
        relative_dir=Plugin.getHome()
        absolute_dir= os.path.join(relative_dir, "config.yml")
        args = f"--config {absolute_dir} --smiles {os.path.abspath(smiles_file)}"
        Plugin.runAIZYNTH("aizynthcli", args, cwd=os.path.abspath(self._getExtraPath()))

        #JSON FILE 
        fnJson = self._getExtraPath("output.json.gz")

        if not os.path.exists(fnJson):
            raise Exception(f"Output JSON not found: {fnJson}")

        with gzip.open(fnJson, 'rb') as f_in:
            with open(self._getExtraPath("output.json"), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        outDir = os.path.abspath(self._getExtraPath())
        jsonPath= os.path.abspath(self._getExtraPath("output.json"))
        args = '{} {}'.format(jsonPath, outDir)
        scriptPath=Plugin.getScriptsPath()
        Plugin.runScript(self, "retrosynthesis_script.py" ,args, env=AIZYNTHFINDER_DIC,cwd=None, scriptDir=scriptPath)
        cleanPath(fnJson)

    def createOutputStep(self):
        file_path = self._getExtraPath("output.json")

        try:
            with open(file_path, 'r') as f:
                data_json = json.load(f)

            trees_info = data_json['data'][0]['trees']

        except FileNotFoundError:
            raise Exception("AiZynthFinder output.json was not found.")

        except (KeyError, IndexError) as e:
            raise Exception("Could not find retrosynthesis routes in output.json: {}".format(e))

        # Collect terminal precursor molecules from one route.
        def getTerminalPrecursors(node):
            precursors = []

            if node.get('type') == 'mol':
                children = node.get('children', [])

                # A molecule without children is a terminal precursor.
                if not children:
                    precursors.append({
                        'smiles': node.get('smiles'),
                        'in_stock': node.get('in_stock', False)
                    })

            # Continue traversing the route tree.
            for child in node.get('children', []):
                precursors.extend(getTerminalPrecursors(child))

            return precursors

        route_outputs = {}
        routes_summary = []

        for route_idx, route in enumerate(trees_info, start=1):

            # Get terminal precursors belonging ONLY to this route.
            precursors = getTerminalPrecursors(route)

            # Remove entries without a SMILES.
            precursors = [
                p for p in precursors if p['smiles']
            ]

            if not precursors:
                print("Route {} has no terminal precursors.".format(route_idx))
                continue

            # Route-level information from AiZynthFinder.
            scores = route.get('scores', {})
            metadata = route.get('metadata', {})

            is_solved = metadata.get('is_solved', False)
            score = scores.get('state score', None)
            number_of_reactions = scores.get('number of reactions', None)

            precursors_in_stock = sum(
                1 for p in precursors if p['in_stock']
            )

            precursors_not_in_stock = len(precursors) - precursors_in_stock

            # Save route information for JSON and CSV summaries.
            route_info = {
                'route': route_idx,
                'solved': is_solved,
                'score': score,
                'reactions': number_of_reactions,
                'precursors': len(precursors),
                'precursors_in_stock': precursors_in_stock,
                'precursors_not_in_stock': precursors_not_in_stock,
                'precursor_smiles': [
                    p['smiles'] for p in precursors
                ]
            }

            routes_summary.append(route_info)

            print("\n================ ROUTE {} ================".format(route_idx))
            print("Solved: {}".format(is_solved))
            print("Score: {}".format(score))
            print("Number of reactions: {}".format(number_of_reactions))
            print("Number of terminal precursors: {}".format(len(precursors)))
            print("Precursors in stock: {}".format(precursors_in_stock))
            print("Precursors not in stock: {}".format(precursors_not_in_stock))

            # Create a separate SetOfSmallMolecules for this route.
            outputRoute = SetOfSmallMolecules().create(
                outputPath=self._getPath(),
                suffix='Route{}'.format(route_idx)
            )

            for precursor_idx, precursor in enumerate(precursors, start=1):

                smile = precursor['smiles']
                in_stock = precursor['in_stock']

                # Generate a unique SDF filename for this route and precursor.
                file_name = "Route{}_Precursor{}".format(
                    route_idx, precursor_idx
                )

                sdf_filename = os.path.abspath(
                    self._getExtraPath(file_name + ".sdf")
                )

                outDir = os.path.abspath(self._getExtraPath())

                args = '--smiles "{}" --output-directory "{}" --filename "{}"'.format(
                    smile, outDir, file_name
                )

                scriptPath = Plugin.getScriptsPath()

                Plugin.runScript(
                    self,
                    "convert_smiles_to_sdf_script.py",
                    args,
                    env=AIZYNTHFINDER_DIC,
                    cwd=None,
                    scriptDir=scriptPath
                )

                # Retrieve molecule name from PubChem.
                cid, link = self.getCIDFromSmiles(smile)
                name = self.getMainNameFromCID(cid)

                if name:
                    name = name.lower()
                else:
                    name = smile

                display_name = self.getSafeMoleculeName(name)

                smallMolecule = SmallMolecule(
                    smallMolFilename=sdf_filename
                )

                smallMolecule.setMolName(display_name)
                outputRoute.append(smallMolecule)

                print(
                    "Precursor {}: {} | {} | SMILES: {}".format(
                        precursor_idx,
                        display_name,
                        link,
                        smile
                    )
                )

            outputRoute.updateMolClass()

            # Register this route as an independent protocol output.
            route_name = 'outputRoute{}'.format(route_idx)
            route_outputs[route_name] = outputRoute

        if not route_outputs:
            raise Exception(
                "No retrosynthesis routes with terminal precursors were found."
            )

        # Register all route outputs.
        self._defineOutputs(**route_outputs)

        # Save route statistics.
        self.saveRoutesSummary(routes_summary)

        print("\nRoute summary saved:")
        print(self._getExtraPath("routes_summary.json"))
        print(self._getPath("routes_summary.csv"))


    
    def getCIDFromSmiles(self, smi):
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/%s/cids/TXT" % smi
        try:
            with urlopen(url) as response:
                cid = response.read().decode('utf-8').split()[0]
        except Exception as e:
            cid = None
        return cid, url
     
    def getMainNameFromCID(self,cid):
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/synonyms/TXT".format(cid)
        try:
            with urlopen(url) as response:
                r = response.read().decode('utf-8')
                synonyms = r.strip().split('\n')
                
                if synonyms:
                    main_name = synonyms[0].strip()
                else:
                    main_name = None
                
        except Exception as e:
            main_name = None
        
        return main_name

    def saveRoutesSummary(self, routes_summary):
        """Save retrosynthesis route statistics in JSON and CSV formats."""

        # JSON: complete information for each route.
        json_file = self._getExtraPath("routes_summary.json")

        with open(json_file, 'w') as f:
            json.dump(routes_summary, f, indent=4)

        # CSV: compact table for easy comparison.
        csv_file = self._getPath("routes_summary.csv")

        fields = [
            'route',
            'solved',
            'score',
            'reactions',
            'precursors',
            'precursors_in_stock',
            'precursors_not_in_stock',
            'precursor_smiles'
        ]
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()

            for route in routes_summary:
                row = {
                    field: route.get(field)
                    for field in fields
                }

                row['precursor_smiles'] = '; '.join(
                    route.get('precursor_smiles', [])
                )

                writer.writerow(row)

    def getSafeMoleculeName(self, name):
        """Replace characters that can cause problems in PyMOL object names."""
        name = re.sub(r'[^a-zA-Z0-9_-]', '_', name)
        name = re.sub(r'_+', '_', name)
        return name.strip('_')

    # --------------------------- INFO functions -----------------------------------


    def _methods(self):
        methods = []
        return methods

    def _validate(self):
         validations = []
         return validations

    def _warnings(self):
         warnings = []
         return warnings

    def _summary(self):
        """Display retrosynthesis statistics in the Scipion Summary tab."""

        summary_file = self._getExtraPath("routes_summary.json")

        if not os.path.exists(summary_file):
            return ["Retrosynthesis route summary is not available yet."]

        with open(summary_file, 'r') as f:
            routes_summary = json.load(f)

        if not routes_summary:
            return ["No retrosynthesis routes were found."]

        solved_routes = sum(
            1 for route in routes_summary if route['solved']
        )

        total_routes = len(routes_summary)

        summary = [
            "Total routes: {}".format(total_routes),
            "Route details:"
        ]

        for route in routes_summary:
            summary.append("")
            summary.append(
                "Route {}: {}".format(
                    route['route'],
                    "Solved" if route['solved'] else "Unsolved"
                )
            )
            summary.append("  Score: {}".format(route['score']))
            summary.append("  Reaction steps: {}".format(route['reactions']))
            summary.append("  Terminal precursors: {}".format(route['precursors']))
            summary.append(
                "  Precursors in stock: {}".format(
                    route['precursors_in_stock']
                )
            )
            summary.append(
                "  Precursors not in stock: {}".format(
                    route['precursors_not_in_stock']
                )
            )

        summary.extend([
            "",
            "Detailed route information:",
            "  routes_summary.json",
            "  routes_summary.csv"
        ])

        return summary

    