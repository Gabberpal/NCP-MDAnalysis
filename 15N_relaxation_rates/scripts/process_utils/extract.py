import csv
import os

from typing import Callable

from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis import AtomGroup


class WriteVectorsToCsv(AnalysisBase):
    def __init__(self,
                 ag: AtomGroup,
                 selector: Callable,
                 filename_provider: Callable,
                 **kwargs):
        """
        :param atomgroup: group of atoms associated with trajectory
        :param selector: selector of atom pairs: atom_1 (origin of vector) and atom_2 (end of vector)
        :param filename_provider: callable object (function) to set output names for files with vector coordinates
        """
        super(WriteVectorsToCsv, self).__init__(ag.universe.trajectory,
                                                **kwargs)
        self._ag = ag

        self.selector = selector
        self.filename_provider = filename_provider
        
        self.files = []
        self.atoms = self.selector(self._ag)

    def _prepare(self):
        os.makedirs(self.filename_provider.output_directory, exist_ok=True)

        for atoms in self.atoms:
            file_name = self.filename_provider(atoms[0], atoms[1])
            file = open(file_name, 'w')
            spamwriter = csv.writer(file)
            spamwriter.writerow(['x', 'y', 'z'])
            self.files.append(file)
        
        self.files_and_atoms = dict.fromkeys(self.files, self.atoms)
            # file = open(file_name, 'w')
            # file.write('x,y,z\n')
            # self.files.append(file_name)
        # called before iteration on the trajectory has begun
        # open csv files to write vectors

    def _single_frame(self):
        for file, atoms in self.files_and_atoms.items():
            spamwriter = csv.writer(file)
            spamwriter.writerow(atoms[0].position - atoms[1].position)
        # write csv vectors for single frame

    def _conclude(self):
        for file in self.files_and_atoms.keys():
            file.close()
        # called once iteration on the trajectory is finished.
        # close csv files
