from typing import Iterable, Callable, List, Tuple, cast
from MDAnalysis import Universe, AtomGroup
from MDAnalysis.core.groups import Atom


def atom_pair_selecetor(atom_name_1: str,
                        atom_name_2: str,
                        resids_of_interest: Iterable[int]) -> Callable:
    """
    :param atom_name_1: name of first atom in pair
    :param atom_name_2: name of second atom in pair
    :param rids_of_interest: set of residue IDs of atoms
    :return: function that allows to select the specified atom pairs
    """

    def select(ag: AtomGroup) -> List[Tuple[Atom, Atom]]:
        """
        :param ag:
        :return:
        """
        return [
            (first_atom, second_atom)
            for resid in resids_of_interest
            for first_atom in
            cast(Iterable[Atom], ag.atoms.select_atoms(f"resid {resid} and name {atom_name_1}"))
            for second_atom in
            cast(Iterable[Atom], ag.atoms.select_atoms(f"resid {resid} and name {atom_name_2}"))
        ]

    return select
