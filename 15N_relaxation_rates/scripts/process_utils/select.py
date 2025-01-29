from typing import Iterable, Callable, List, Tuple, cast
from MDAnalysis import Universe, AtomGroup
from MDAnalysis.core.groups import Atom
from MDAnalysis.analysis import dssp


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


def get_sec_str_pattern(reference: Universe,
                        cnain_ids: Iterable[str]
                        ) -> str:
    """
    :param reference: srtucture to assign secondary-structured elements
    :param chainIDs: list of chain names: ["A"] or ["A", "B"]
    :return: string pattern for selection of atoms from secondary structured regions in specified chains
    """
    sec_str_patterns = []
    for chain in cnain_ids:

        dssp_results = dssp.DSSP(reference.select_atoms(f"chainID {chain}")).run().results
        ss_indexes = dssp_results.dssp_ndarray[0][:, 1:].sum(axis=1).astype(bool)
        ss_residues = dssp_results.resids[ss_indexes]

        if len(ss_residues) > 0:
            ss_selection = f"(chainID {chain} and resid {' '.join(ss_residues.astype(str))})"
            sec_str_patterns.append(ss_selection)

    return f"({' or '.join(sec_str_patterns)})" if sec_str_patterns else ""
