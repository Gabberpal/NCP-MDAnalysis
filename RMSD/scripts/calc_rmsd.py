import argparse
import os
import pandas as pd
from glob import glob

import MDAnalysis as mda
import MDAnalysis.transformations as trans
from MDAnalysis import Universe
from MDAnalysis.analysis import rms

from process_utils.select import get_sec_str_pattern

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
    trj_list = glob(os.path.join(args.path_to_trajectory, "*nc"))
    trj_list.sort()
    u = Universe(args.path_to_trajectory_reference,
                 trj_list[args.trajectory_start:args.trajectory_length],
                 in_memory=True,
                 in_memory_step=args.frames_per_trajectory_file
                 )
    ref_trj = mda.Universe(args.path_to_trajectory_reference)
    chainids = []
    for segment in u.segments:
        chainids.extend([segment.segid] * segment.atoms.n_atoms)
    u.add_TopologyAttr("chainIDs", chainids)
    ref_trj.add_TopologyAttr("chainIDs", chainids)

    # set xray reference to calc RMSD
    xray_reference = mda.Universe(args.path_to_xray_reference)

    # set pattern to select CA atoms from secondary structure
    chainids = []
    for segment in xray_reference.segments:
        chainids.extend([segment.segid] * segment.atoms.n_atoms)
    xray_reference.add_TopologyAttr("chainIDs", chainids)

    protein_chains = ["A", "B", "C", "D", "E", "F", "G", "H"]
    selection_sec_str = get_sec_str_pattern(reference=xray_reference,
                                            cnain_ids=protein_chains)
    selection_sec_str_ca = f"name CA and {selection_sec_str}"

    # set pattern to select inner and outer DNA turns
    dna_inner_turn = ' '.join(f'{i}' for i in range(-38, 38 + 1))
    dna_outer_turn = ' '.join(f'{i}' for i in range(-72, -39 + 1)) + f' ' \
                     + ' '.join(f'{i}' for i in range(39, 72 + 1))

    inner_dna_seceletion = f"(name N1 or name N9) and (chainID I or chainID J) and (resid {dna_inner_turn})"
    outer_dna_seceletion = f"(name N1 or name N9) and (chainID I or chainID J) and (resid {dna_outer_turn})"

    dna = f"({inner_dna_seceletion}) or ({outer_dna_seceletion})"
    all = f"({dna}) or ({selection_sec_str})"

    # transform trajectory
    atoms = u.atoms
    transforms = [
        trans.NoJump(),
        trans.center_in_box(atoms),
        trans.wrap(atoms, compound="segments"),
        trans.fit_rot_trans(atoms, ref_trj)
    ]
    u.trajectory.add_transformations(*transforms)

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