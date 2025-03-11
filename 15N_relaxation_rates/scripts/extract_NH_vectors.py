import argparse
import os

from glob import glob

from MDAnalysis import Universe
from MDAnalysis.core.groups import Atom

from process_utils.select import atom_pair_selecetor, get_sec_str_ca_pattern
from process_utils.analysis import ExtractVectors, AnalyzerWrapper
from process_utils.transform import AssembleQuaternaryStructure, TransformWrapper, fit_rot_trans_by_pattern
from process_utils.batch_process import BatchLoader, BatchCsvWriter, BatchAnalyzer


def vectorname_provider(atom_1: Atom, atom_2: Atom):
    return f"{atom_1.chainID}_{atom_1.resid:02d}_{atom_1.resname}_{atom_1.name}{atom_2.name}"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract NH vectors')
    parser.add_argument('--path-to-trajectory', required=True)
    parser.add_argument('--path-to-trajectory-reference', required=True)
    parser.add_argument('--chain-name', required=True)
    parser.add_argument('--residue-of-interest', required=True)
    parser.add_argument('--trajectory-start', default=0, type=int)
    parser.add_argument('--trajectory-length', required=True, type=int)
    parser.add_argument('--trajectory-stride', default=1, type=int)
    parser.add_argument('--batch-size', default=1, type=int)
    parser.add_argument('--dt-ns', type=float, default=0.01)
    parser.add_argument('--output-directory', default=".")
    parser.add_argument('--dna-chains', default=["I", "J"], type=list)
    parser.add_argument('--protein-chains', default=["A", "B", "C", "D", "E", "F", "G", "H"], type=list)
    args = parser.parse_args()

    #  load trajectory reference
    trj_reference = Universe(args.path_to_trajectory_reference, topology_format="PDB")

    # set path to trajectory files
    nc_files = glob(os.path.join(args.path_to_trajectory, '*.nc'))
    nc_files.sort()

    # set trajectory transforms
    transforms = [
        # 1. assemble the nucleosome particle that may appear divided at the boundaries
        # due to periodic boundary condition in the MD simulation.
        TransformWrapper(transform=AssembleQuaternaryStructure,
                         reference=trj_reference,
                         cnain_ids=args.protein_chains + args.dna_chains,
                         atom_selector="name CA N1 N9"),
        # 2. overlay all MD frames by superimposing them onto the reference
        # via the secondary-structure Cα atoms from the histone core.
        TransformWrapper(transform=fit_rot_trans_by_pattern,
                         reference=trj_reference,
                         pattern=get_sec_str_ca_pattern(reference=trj_reference,
                                                        cnain_ids=args.protein_chains
                                                        ))
    ]

    # set loader for trajectory to batch processing due to file open limits
    batchloader = BatchLoader(reference=trj_reference,
                              trj_list=nc_files[args.trajectory_start:args.trajectory_length],
                              trajectory_stride=args.trajectory_stride,
                              batch_size=args.batch_size,
                              dt_ns=args.dt_ns,
                              transforms=transforms,
                              pattern=f"chainID {args.chain_name}"
                              )

    # set selector to extract N-H vectors from trajectory
    first_rid, last_rid = args.residue_of_interest.split("-")
    resids_of_interest = set(list(range(int(first_rid), int(last_rid) + 1)))
    selector = atom_pair_selecetor(atom_name_1="N",
                                   atom_name_2="H",
                                   resids_of_interest=resids_of_interest)

    # set trajectory analyzer to extract coordinates of NH vectors
    analyzer = AnalyzerWrapper(ExtractVectors,
                               selector=selector,
                               vectorname_provider=vectorname_provider,
                               )

    # set writer to save the coordinates of NH vectors
    writer = BatchCsvWriter(output_directory=args.output_directory,
                            header=["time_ns", "x", "y", "z"]
                            )

    # process batches of trajectory and save the results
    batchanalyzer = BatchAnalyzer(batchloader=batchloader,
                                  analyzer=analyzer,
                                  writer=writer)
    batchanalyzer.analyse()