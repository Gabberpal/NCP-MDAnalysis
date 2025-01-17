import argparse
import MDAnalysis as mda

from MDAnalysis import Universe
from process_utils.select import get_sec_str_pattern

import pandas as pd
import os
from MDAnalysis.analysis import rms


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculation RMSD')
    parser.add_argument('--path-to-trajectory', required=True)
    parser.add_argument('--path-to-trajectory-reference', required=True)
    parser.add_argument('--path-to-xray-reference', required=True)
    parser.add_argument('--filetype', choices=["dat", "nc", "xtc"], required=True)
    parser.add_argument('--pattern', default="run%05d")
    parser.add_argument('--trajectory-start', default=1, type=int)
    parser.add_argument('--trajectory-length', required=True, type=int)
    parser.add_argument('--trajectory-stride', default=1, type=int)
    parser.add_argument('--frames-per-trajectory-file', type=int, default=100)
    parser.add_argument('--dt-ns', type=float, default=0.01)
    parser.add_argument('--output-directory', default=".")
    args = parser.parse_args()

    #  load trajectory
    u: Universe = ...  # TODO: Implement reading trajectory

    # set xray reference to calc RMSD
    xray_reference = mda.Universe(args.path_to_xray_reference)

    # set pattern to select CA atoms from secondary structure
    protein_chains = "A", "B", "C", "D", "E", "F", "G", "H"
    # TODO: get_sec_str_pattern
    selection_sec_str = get_sec_str_pattern(reference=xray_reference,
                                            cnain_ids=protein_chains)
    selection_sec_str_ca = f"name CA and {selection_sec_str}"

    # set pattern to select inner and outer DNA turns
    # TODO:
    ...

    # transform trajectory
    # TODO:
    ...

    # calculate and save RMSD data
    R = rms.RMSD(atomgroup=u,
                 reference=ref_trj,
                 select=selection_sec_str,
                 groupselections=[inner_dna_seceletion, outer_dna_seceletion, dna, all],
                 ref_frame=0)
    R.run()


    rmsd_df = pd.DataFrame(R.results.rmsd[:, 1:],
                           columns=['time_ns', 'rmsd_protein', 'rmsd_dna_inner', 'rmsd_dna_outer', 'rmsd_dna',
                                    'rmsd_all'])
    rmsd_df["time_ns"] /= 1000


    # saving to a CSV file
    os.makedirs(args.output_directory, exist_ok=True)
    rmsd_df.to_csv(os.path.join(args.output_directory, "rmsd.csv"), index=False)