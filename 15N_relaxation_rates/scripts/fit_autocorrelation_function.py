from process_utils.save_utils import fit_and_save_acorr_func
import argparse
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='fit acorr')
    parser.add_argument('--path-to-acorr-csv', required=True)
    parser.add_argument('--limit_ns', default=None, type=int)
    parser.add_argument('--lag-spacing', default=None)
    parser.add_argument('--n-lag-points', default=1000, type=int)
    parser.add_argument('--path-to-reference-pdb', type=str)
    parser.add_argument('--chain-name', type=str)
    parser.add_argument('--fast-timescale-ns', default=0.0001, type=float)
    parser.add_argument('--tumbling-time-ns', type=float)
    parser.add_argument('--number_of-exponent', default=6, type=int)
    parser.add_argument('--output-directory', default="./")
    args = parser.parse_args()

    n_exp = args.number_of_exponent
    fast_timescale = args.fast_timescale_ns
    tumbling_time = args.tumbling_time_ns

    bounds = [[[0, fast_timescale] * n_exp,
               [1.0, tumbling_time] * n_exp]]

    p0 = np.zeros(n_exp * 2)
    p0[::2] = np.array([1.0 / n_exp] * n_exp)
    p0[1::2] = np.logspace(np.log10(fast_timescale), np.log10(round(tumbling_time, 0)), n_exp)

    fit_and_save_acorr_func(args.path_to_acorr_csv,
                            bounds=bounds,
                            p0=p0[1:],
                            limit_ns=args.limit_ns,
                            lag_spacing=args.lag_spacing,
                            n_lag_points=args.n_lag_points,
                            output_directory=args.output_directory
                            )