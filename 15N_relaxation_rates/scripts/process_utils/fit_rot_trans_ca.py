from MDAnalysis.transformations import fit_rot_trans
from MDAnalysis.analysis import align


class fit_rot_trans_ca(fit_rot_trans):

    def __init__(self, ag, reference, plane=None, weights=None, max_threads=1, parallelizable=True):
        super().__init__(ag, reference, plane=plane, weights=weights, max_threads=max_threads,
                         parallelizable=parallelizable)

        self.mobile_all_atoms_com = self.mobile.atoms.center(self.weights)

        ag_ca = ag.select_atoms("name CA")
        ref_ca = reference.select_atoms("name CA")
        self.ref_com_all = self.ref.center(self.weights)

        self.ref, self.mobile = align.get_matching_atoms(ref_ca.atoms,
                                                         ag_ca.atoms)

        self.weights = align.get_weights(self.ref.atoms,
                                         weights=self.weights)

        self.ref_com = self.ref.center(self.weights)
        self.ref_coordinates = self.ref.atoms.positions - self.ref_com
