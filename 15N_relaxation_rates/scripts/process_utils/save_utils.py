import os
import pandas as pd
import numpy as np 

from glob import glob
from tqdm import tqdm

from typing import Callable, Iterable
from process_utils.calc import calc_acorr_order_2
from process_utils.fit import repeated_fit_auto_correlation


def calc_and_save_acorr(path_to_vector_csv_files: Iterable[str],
                        dt_ns: float,
                        acorr_func_limit: int = -1,
                        thumbling_time: float = None,
                        acorr_func: Callable = calc_acorr_order_2,
                        out_dir: str = "."):
    """
    This function calculates the autocorrelation function for a vector array stored in CSV files and saves the results to an output directory.

    :param path_to_vector_csv_files: Path to a single CSV file or a list of CSV file paths containing the 3-dimensional vector arrays [N x 3].
    :param dt_ns: Time step (in nanoseconds) between consecutive data points in the vector array
    :param acorr_func_limit: The maximum lag (in terms of data points) for which the autocorrelation is stored.
                             If set to -1 (default), the autocorrelation is stored for all possible lags.
    :param thumbling_time: if is not None the overall tumbling is reintroduced by multiplying the result
                           by exp(-τ/τ_rot), where τ_rot is rotational correlation time of the protein
                           Default is None
    :param acorr_func: The function used to calculate the autocorrelation.
                       This function should accept a 3D numpy array and return the autocorrelation values.
                       Default is calc_acorr_order_2.
    :param out_dir: Directory where the autocorrelation results will be saved. Default is the current directory
    :return:
    """
    index = None
    for path_to_vector_file in tqdm(sorted(path_to_vector_csv_files)):

        out_name = os.path.basename(path_to_vector_file).split("_")[0]
        vectors = pd.read_csv(path_to_vector_file).values
        acorr = acorr_func(vectors)[:acorr_func_limit]

        if index is None:
            index = len(acorr)
            time_ns = np.linspace(0, index * dt_ns, index, endpoint=False)

        if (thumbling_time is not None) and (thumbling_time > 0):
            acorr *= np.exp(-time_ns / thumbling_time)

        os.makedirs(out_dir, exist_ok=True)
        pd.DataFrame({
            "time_ns": time_ns,
            "acorr": acorr
        }).to_csv(os.path.join(out_dir, out_name), index=False)


def fit_and_save_acorr_func(path_to_acorr_files,
                            bounds,
                            p0=None,
                            lag_spacing="log",
                            n_lag_points=None,
                            output_directory="./",
                            limit_ns=None):
    path_to_ccr_csv_files = sorted(glob(os.path.join(path_to_acorr_files, "*.csv")))
    rname_list = []
    for file in path_to_ccr_csv_files:
        aa_name = file.split('-')[2]
        rname_list.append(aa_name)
    for bound in bounds:
        tau_table = pd.DataFrame()
        for acorr_corr_file in tqdm(path_to_ccr_csv_files, desc=output_directory):
            df = pd.read_csv(acorr_corr_file)
            if limit_ns:
                df = df[df["time_ns"] <= limit_ns]

            time_ns, acorr = df["time_ns"].values, df["acorr"].values

            if lag_spacing == "log":
                lag_index = np.unique(
                    np.logspace(0, int(np.log10(time_ns.size)), n_lag_points, endpoint=False).astype(int))
                acorr = np.take(acorr, lag_index)
                time_ns = np.take(time_ns, lag_index)

            popt = repeated_fit_auto_correlation(acorr, time_ns, bound, p0)
            name = os.path.splitext(os.path.basename(acorr_corr_file))[0]
            amplitudes = popt[::2]
            taus = popt[1::2]
            order = (len(bound[0]) + 1) // 2

            rid = int(name.split("-")[1])

            popt_dict = {
                'rId': rid,
                'rName': rname_list[rid - 2],
                'limit_ns': limit_ns
            }

            popt_dict.update(
                {("exp-%d-a%d" % (order, i + 1)): a for i, a in enumerate(amplitudes)}
            )
            popt_dict.update(
                {("exp-%d-tau%d" % (order, i + 1)): tau for i, tau in enumerate(taus)}
            )

            tau_table = pd.concat([tau_table, pd.DataFrame(popt_dict, index=[0])])

        tau_table.sort_values(by=['rId'], inplace=True)
        os.makedirs(output_directory, exist_ok=True)
        tau_table.to_csv(os.path.join(output_directory, 'tau_%d_exp.csv' % order), index=False)