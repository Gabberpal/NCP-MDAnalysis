from typing import Iterable
from MDAnalysis import Universe


def get_sec_str_pattern(reference: Universe,
                        cnain_ids: Iterable[str]
                        ) -> str:
    """
    :param reference: srtucture to assign secondary-structured elements
    :param chainIDs: list of chain names: ["A"] or ["A", "B"]
    :return: string pattern for selection of atoms from secondary structured regions in specified chains
    """
    pass
