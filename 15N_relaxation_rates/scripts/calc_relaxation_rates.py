import argparse

from process_utils.calc_relaxation_rate import get_relaxition_rate

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calc DD-DD relaxation rate')
    parser.add_argument('--path-to-fit-dir', required=True, )
    parser.add_argument('--nmr-freq', required=True, type=float)
    parser.add_argument('--output-directory', default="./")

    args = parser.parse_args()

    for rate in ["R1", "R2"]:
        get_relaxition_rate(path_to_fit=args.path_to_fit_dir,
                            nmr_freq=args.nmr_freq, rate=rate,
                            output_directory=args.output_directory)
