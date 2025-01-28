from process_utils.save_utils import plot_and_save_acorr_with_fit
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='fit acorr')
    parser.add_argument('--path-to-acorr-csv', required=True)
    parser.add_argument('--path-to-fit-csv', required=True)
    parser.add_argument('--path-to-reference-pdb', required=True)
    parser.add_argument('--output-directory', default="./")
    args = parser.parse_args()

    plot_and_save_acorr_with_fit(path_to_fit_csv=args.path_to_fit_csv,
                                 path_to_acorr_csv=args.path_to_acorr_csv,
                                 output_directory=args.output_directory
                                 )
