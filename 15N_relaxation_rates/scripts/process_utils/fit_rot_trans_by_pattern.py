from MDAnalysis.transformations import fit_rot_trans
from MDAnalysis.analysis import align


class fit_rot_trans_by_pattern(fit_rot_trans):

    def __init__(self, ag, reference, pattern=None, *args, **kwargs):
        super().__init__(ag, reference, *args, **kwargs)

        selection_pattern = pattern if pattern is not None else 'all'

        self.mobile_all_atoms_com = self.mobile.atoms.center(self.weights)

        ag_ca = ag.select_atoms(selection_pattern)
        ref_ca = reference.select_atoms(selection_pattern)
        self.ref_com_all = self.ref.center(self.weights)

        self.ref, self.mobile = align.get_matching_atoms(ref_ca.atoms,
                                                         ag_ca.atoms)

        self.weights = align.get_weights(self.ref.atoms,
                                         weights=self.weights)

        self.ref_com = self.ref.center(self.weights)
        self.ref_coordinates = self.ref.atoms.positions - self.ref_com

