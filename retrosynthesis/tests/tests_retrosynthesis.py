
#### VERÓNICA GAMO ######

# Scipion chem imports
from retrosynthesis.protocols import ProtChemAiZynthFinder
from pwchem.utils import assertHandle
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pwchem.protocols import ProtChemImportSmallMolecules


class TestImportBase(BaseTest):
	@classmethod
	def setUpClass(cls):
		cls.dsLig = DataSet.getDataSet("smallMolecules")
		setupTestProject(cls)

		cls._runImportSmallMols()
		cls._waitOutput(cls.protImportSmallMols, 'outputSmallMolecules', sleepTime=5)

	@classmethod
	def _runImportSmallMols(cls):
		cls.protImportSmallMols = cls.newProtocol(
			ProtChemImportSmallMolecules,
			filesPath=cls.dsLig.getFile('sdf'))
		cls.proj.launchProtocol(cls.protImportSmallMols, wait=False)

class TestRetrosynthesis(TestImportBase):
	@classmethod
	def _runAizynthfinder(cls, inProt, mode=0):
		protaizynth = cls.newProtocol(
			ProtChemAiZynthFinder,
			inputLigand='SmallMolecule (2000 molecule)'
		)
		protaizynth.inputSet.set(inProt)
		protaizynth.inputSet.setExtended('outputSmallMolecules')

		cls.proj.launchProtocol(protaizynth, wait=False)
		return protaizynth

	def test(self):
		p= self._runAizynthfinder(inProt=self.protImportSmallMols)
		self._waitOutput(p, 'outputSmallMolecules', sleepTime=10)
		assertHandle(self.assertIsNotNone, getattr(p, 'outputSmallMolecules', None), cwd=p.getWorkingDir())


