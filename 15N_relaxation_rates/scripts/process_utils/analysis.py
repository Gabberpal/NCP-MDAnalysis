import numpy as np
from typing import Callable

from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis import AtomGroup


class AnalyzerWrapper:
    def __init__(self, analysis, *args, **kwargs):
        self.analysis = analysis

        self.__dict__.update(kwargs)
        self.attr_names = list(kwargs.keys())

        for arg in args:
            setattr(self, arg, arg)
            self.attr_names.append(arg)

    def __call__(self, ag):
        attrs = [getattr(self, arg) for arg in self.attr_names]
        return self.analysis(ag, *attrs)


class ExtractVectors(AnalysisBase):
    def __init__(self,
                 ag: AtomGroup,
                 selector: Callable,
                 vectorname_provider: Callable,
                 **kwargs):
        """
        :param atomgroup: group of atoms associated with trajectory
        :param selector: selector of atom pairs: atom_1 (origin of vector) and atom_2 (end of vector)
        :param filename_provider: callable object (function) to set output names for files with vector coordinates
        """
        super(ExtractVectors, self).__init__(ag.universe.trajectory,
                                             **kwargs)
        self._ag = ag

        self.selector = selector
        self.vectorname_provider = vectorname_provider

        self.atom_pairs = self.selector(self._ag)
        self.vector_names = []

    def _prepare(self):
        for atom1, atom2 in self.atom_pairs:
            atom_pair_name = self.vectorname_provider(atom1, atom2)
            self.results.update({atom_pair_name: np.zeros((self.n_frames, 5))})
            self.vector_names.append(atom_pair_name)

    def _single_frame(self):

        for vector_name, (atom1, atom2) in dict(zip(self.vector_names, self.atom_pairs)).items():
            self.results[vector_name][self._frame_index, :2] = self._ts.frame, self._trajectory.time
            self.results[vector_name][self._frame_index, 2:] = atom1.position - atom2.position
