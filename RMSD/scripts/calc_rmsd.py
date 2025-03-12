import argparse
import os

from itertools import chain
from glob import glob

from MDAnalysis import Universe

from process_utils.select import get_sec_str_ca_pattern
from MDAnalysis.analysis.rms import RMSD
from process_utils.transform import AssembleQuaternaryStructure, TransformWrapper
from process_utils.batch_process import BatchLoader, BatchCsvWriter, BatchAnalyzer


class AnalyzerWrapper:
    def __init__(self, analysis, *args, **kwargs):
        self.analysis = analysis

        self.__dict__.update(kwargs)
        self.attr_names = list(kwargs.keys())

        for arg in args:
            setattr(self, arg, arg)
            self.attr_names.append(arg)

    def __call__(self, ag):
        attrs = [getattr(self, arg) for arg in self.attr_names]
        return self.analysis(ag, *attrs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculation RMSD')
    parser.add_argument('--path-to-trajectory', required=True)
    parser.add_argument('--path-to-trajectory-reference', required=True)
    parser.add_argument('--path-to-xray-reference', required=True)
    parser.add_argument('--trajectory-start', default=0, type=int)
    parser.add_argument('--trajectory-length', required=True, type=int)
    parser.add_argument('--trajectory-stride', default=100, type=int)
    parser.add_argument('--batch-size', default=1, type=int)
    parser.add_argument('--dt-ns', type=float, default=0.01)
    parser.add_argument('--output-directory', default=".")
    parser.add_argument('--dna-chains', default=["I", "J"], type=list)
    parser.add_argument('--protein-chains', default=["A", "B", "C", "D", "E", "F", "G", "H"], type=list)
    args = parser.parse_args()

    # set xray reference to calc RMSD
    xray_reference = Universe(args.path_to_xray_reference, topology_format="PDB")

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
    ]

    # set loader for trajectory to batch processing due to file open limits
    batchloader = BatchLoader(reference=trj_reference,
                              trj_list=nc_files[args.trajectory_start:args.trajectory_length],
                              trajectory_stride=args.trajectory_stride,
                              batch_size=args.batch_size,
                              dt_ns=args.dt_ns,
                              transforms=transforms,
                              pattern="name CA N1 N9"
                              )

    # set calculate RMSD using different sets of atoms:
    # (1) Cα atoms that have been used to overlay the frames
    # (2) N1 and N9 atoms from the inner turn of nDNA, nucleotides from -38 to 38
    # (3) N1 and N9 atoms from the outer turn of nDNA, nucleotides from -72 to -39 and from 39 to 72
    # (4) N1 and N9 atoms from nDNA nucleobases
    # (5) sets (1) and (4) combined
    sec_str_ca = get_sec_str_ca_pattern(xray_reference, cnain_ids=args.protein_chains)

    dna_inner_turn = " ".join(f'{i}' for i in range(-38, 38 + 1))
    dna_outer_turn = ' '.join(f'{i}' for i in chain(range(-72, -39 + 1), range(39, 72 + 1)))

    inner_dna_pattern = f"(chainID {' '.join(args.dna_chains)}) and (resid {dna_inner_turn}) and (name N1 N9)"
    outer_dna_pattern = f"(chainID {' '.join(args.dna_chains)})  and (resid {dna_outer_turn}) and (name N1 N9)"

    dna_pattern = f"({inner_dna_pattern}) or ({outer_dna_pattern})"
    sec_str_ca_and_dna_pattern = f"({dna_pattern}) or ({sec_str_ca})"

    # set trajectory analyzer to extract coordinates of NH vectors
    analyzer = AnalyzerWrapper(RMSD,
                               reference=xray_reference,
                               select=sec_str_ca,
                               groupselections=[dna_pattern, sec_str_ca_and_dna_pattern,
                                                inner_dna_pattern, outer_dna_pattern]
                               )

    # set writer to save the coordinates of NH vectors
    writer = BatchCsvWriter(output_directory=args.output_directory,
                            header=["time_ns", "rmsd_protein", "rmsd_dna", "rmsd_all",
                                    "rmsd_dna_inner", "rmsd_dna_outer"]
                            )

    # process batches of trajectory and save the results
    batchanalyzer = BatchAnalyzer(batchloader=batchloader,
                                  analyzer=analyzer,
                                  writer=writer)
    batchanalyzer.analyse()
