import argparse

from MDAnalysis import Universe
from MDAnalysis.core.groups import Atom

from process_utils.select import atom_pair_selecetor
from process_utils.extract import WriteVectorsToCsv


class OutputFilenameFormatter:
    def __init__(self, output_directory):
        self.output_directory = output_directory

    def __call__(self, atom_1: Atom, atom_2: Atom):
        return f"{self.output_directory}/" \
               f"{atom_1.chainID}-{atom_1.resid:02d}-{atom_1.resname}" \
               f"-NH.csv"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract NH vectors')
    parser.add_argument('--path-to-trajectory', required=True)
    parser.add_argument('--path-to-trajectory-reference', required=True)
    parser.add_argument('--chain-name', required=True)
    parser.add_argument('--residue-of-interest', required=True)
    parser.add_argument('--filetype', choices=["nc", "xtc"], required=True)
    parser.add_argument('--pattern', default="run%05d")
    parser.add_argument('--trajectory-start', default=1, type=int)
    parser.add_argument('--trajectory-length', required=True, type=int)
    parser.add_argument('--output-directory', default=".")
    args = parser.parse_args()

    #  load trajectory
    traj = sorted(glob(os.path.join(args.path_to_trajectory, '*.nc'))) 
    ref = Universe(args.path_to_trajectory_reference, format='PDB')
    u = Universe(args.path_to_trajectory_reference, traj[:], in_memory=True, topology_format='PDB')

    # transform trajectory:
    transforms = [
        trans.NoJump(),
        trans.center_in_box(u.atoms),
        trans.wrap(u.atoms, compund='segments'),
        fit_rot_trans_ca(u.atoms, ref)
    ]
    u.trajectory.add_transformations(*transforms)

    # write coordinates of N-H vectors for residues of interest
    first_rid, last_rid = args.residue_of_interest.split("-")
    resids_of_interest = set(list(range(int(first_rid), int(last_rid) + 1)))

    # TODO: Implement WriteVectorsToCsv
    WriteVectorsToCsv(ag=u.select_atoms(f"chainID {args.chain_name}"),
                      selector=atom_pair_selecetor(atom_name_1="N",
                                                   atom_name_2="H",
                                                   resids_of_interest=resids_of_interest),
                      filename_provider=OutputFilenameFormatter(output_directory="."),
                      ).run(verbose=True, progressbar_kwargs={"desc": "extract N-H vectors"})
