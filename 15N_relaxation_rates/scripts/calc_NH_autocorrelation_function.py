import os
import argparse
from glob import glob
from process_utils.save_utils import calc_and_save_acorr

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculation N-H autocorrelation function')
    parser.add_argument('--path-to-vectors-csv-dir', required=True)
    parser.add_argument('--dt-ns', default=0.001, type=float)
    parser.add_argument('--tumbling-time-ns', type=float)
    parser.add_argument('--acorr-func-limit', default=-1, type=int)
    parser.add_argument('--output-directory', default=".")
    args = parser.parse_args()

    vectors_csv = glob(os.path.join(args.path_to_vectors_csv_dir, "*csv"))

    os.makedirs(args.output_directory, exist_ok=True)

    calc_and_save_acorr(vectors_csv,
                        acorr_func_limit=args.acorr_func_limit,
                        dt_ns=args.dt_ns,
                        thumbling_time=args.tumbling_time_ns,
                        out_dir=args.output_directory)
