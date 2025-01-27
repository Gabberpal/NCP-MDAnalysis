from typing import Iterable
from MDAnalysis import Universe
from MDAnalysis.analysis import dssp


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

    return " or ".join(sec_str_patterns) if sec_str_patterns else ""
