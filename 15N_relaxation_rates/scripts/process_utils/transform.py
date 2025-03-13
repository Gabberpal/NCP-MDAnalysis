import numpy as np

from MDAnalysis.analysis import align
from MDAnalysis.transformations import fit_rot_trans
from MDAnalysis.transformations.base import TransformationBase


class TransformWrapper:
    def __init__(self, transform, *args, **kwargs):
        self.transfom = transform

        self.__dict__.update(kwargs)
        self.attr_names = list(kwargs.keys())

        for arg in args:
            setattr(self, arg, arg)
            self.attr_names.append(arg)

    def __call__(self, ag):
        attrs = [getattr(self, arg) for arg in self.attr_names]
        return self.transfom(ag, *attrs)


def metric_tensor(v1, v2, v3):
    basis = np.array([v1, v2, v3])
    g_ij = np.zeros((3, 3))

    for i in range(3):
        for j in range(3):
            g_ij[i, j] = np.dot(basis[i], basis[j])
    return g_ij


class Alignment:

    def __init__(self, mobile, ref, weights=None):
        self.mobile = mobile
        self.ref = ref
        self.weights = weights

    @property
    def matrix3d(self):
        mobile0 = self.mobile.positions - self.mobile.atoms.center_of_mass()
        ref0 = self.ref.positions - self.ref.atoms.center_of_mass()
        R, _ = align.rotation_matrix(mobile0, ref0, weights=self.weights)
        return R

    @property
    def vector3d(self):
        return -np.dot(self.mobile.center_of_mass(), self.matrix3d.T) + self.ref.center_of_mass()


class AssembleQuaternaryStructure(TransformationBase):

    def __init__(self,
                 ag,
                 reference,
                 cnain_ids,
                 atom_selector,
                 max_threads=1,
                 parallelizable=True
                 ):

        super().__init__(
            max_threads=max_threads, parallelizable=parallelizable
        )

        self.ag = ag
        self.reference = reference
        self.atom_selector = atom_selector

        self.reference_mols = [
            self.reference.select_atoms(f"chainID {chain_id}") for chain_id in cnain_ids
        ]

        self._reference_mean_coords = [
            mol.atoms.select_atoms(atom_selector).center_of_mass()
            for mol in self.reference_mols
        ]

        self._mols = [
            self.ag.select_atoms(f"chainID {chain_id}") for chain_id in cnain_ids
        ]

    def _transform(self, ts):
        # first molecule in assembly is aligned by convention
        alignment = Alignment(self.reference_mols[0].select_atoms(self.atom_selector),
                              self._mols[0].select_atoms(self.atom_selector))

        # calculate reference coordinates for current frame
        ref_mean_coords = [np.dot(coords, alignment.matrix3d.T) + alignment.vector3d
                           for coords in self._reference_mean_coords]

        # shift rest of molecules to match reference coordinates as close as possible
        for mol_n in range(1, len(ref_mean_coords)):
            ref_point = ref_mean_coords[mol_n]
            mol = self._mols[mol_n]
            mol_point = self._mols[mol_n].select_atoms(self.atom_selector).center_of_mass()
            v1, v2, v3 = ts.triclinic_dimensions

            # find the closest image
            g_ij = metric_tensor(v1, v2, v3)
            delta = ref_point - mol_point

            approx = np.linalg.inv(g_ij) @ np.array([delta.dot(v1), delta.dot(v2), delta.dot(v3)])
            i, j, k = map(round, approx)

            if sum([i, j, k]):
                shift = v1 * i + v2 * j + v3 * k
                ts.positions[mol.ix] += shift
        return ts


class fit_rot_trans_by_pattern(fit_rot_trans):

    def __init__(self, ag, reference, pattern='all', *args, **kwargs):
        super().__init__(ag, reference, *args, **kwargs)

        self.mobile_all_atoms_com = self.mobile.atoms.center(self.weights)

        ag_ca = ag.select_atoms(pattern)
        ref_ca = reference.select_atoms(pattern)
        self.ref_com_all = self.ref.center(self.weights)

        self.ref, self.mobile = align.get_matching_atoms(ref_ca.atoms,
                                                         ag_ca.atoms)

        self.weights = align.get_weights(self.ref.atoms,
                                         weights=self.weights)

        self.ref_com = self.ref.center(self.weights)
        self.ref_coordinates = self.ref.atoms.positions - self.ref_com
