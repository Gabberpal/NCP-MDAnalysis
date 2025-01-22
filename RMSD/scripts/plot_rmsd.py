from process_utils.save_utils import plot_and_save_rmsd
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='plot RMSD')
    parser.add_argument('--path-to-rmsd_dir', required=True)
    parser.add_argument('--output-directory', default="./")
    args = parser.parse_args()

    plot_and_save_rmsd(path_to_rmsd_dir=args.path_to_rmsd_dir,
                       output_directory=args.output_directory
                       )
