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

    def _prepare(self):
        # called before iteration on the trajectory has begun
        # open csv files to write vectors
        pass

    def _single_frame(self):
        # write csv vectors for single frame
        pass

    def _conclude(self):
        # called once iteration on the trajectory is finished.
        # close csv files
        pass
