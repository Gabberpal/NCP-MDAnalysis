from typing import Iterable
from MDAnalysis import Universe
from MDAnalysis.analysis import dssp


def get_sec_str_pattern(reference: Universe,
                        cnain_ids: Iterable[str]
                        ) -> str:

    select_ca_parts = []
    for protein_chain in cnain_ids:

        results = dssp.DSSP(reference.select_atoms(f"chainID {protein_chain}")).run().results
        ss_indexes = results.dssp_ndarray[0][:, 1:].sum(axis=1).astype(bool)
        ss_residues = results.resids[ss_indexes]

        if len(ss_residues) > 0:
            resid_string = " ".join(ss_residues.astype(str))
            ss_selection = f"(chainID {protein_chain} and resid {resid_string})"
            select_ca_parts.append(ss_selection)
    return " or ".join(select_ca_parts) if select_ca_parts else ""