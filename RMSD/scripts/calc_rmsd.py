import argparse
import MDAnalysis as mda

from MDAnalysis import Universe
from process_utils.select import get_sec_str_pattern

import os
import MDAnalysis.transformations as trans
from glob import glob

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
    traj = Universe(args.path_to_reference_pdb,
                 trj_list[args.trajectory_start:args.trajectory_length + 1 ],
                 in_memory=True,
                 in_memory_step=args.frames_per_trajectory_file,
                 topology_format = "PDB"
                 )

    # set xray reference to calc RMSD
    path_to_ref = args.path_to_xray_reference_pdb
    ref_trj = mda.Universe(path_to_ref, topology_format="PDB")

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
    atoms = u.atoms
    transforms = [
        trans.NoJump(),
        trans.center_in_box(atoms),
        trans.wrap(atoms, compound="segments"),
        trans.fit_rot_trans(atoms, ref_trj)
    ]
    u.trajectory.add_transformations(*transforms)

    # calculate and save RMSD data
    # TODO:
    ...
