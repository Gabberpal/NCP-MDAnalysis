import os
import pandas as pd
import numpy as np 

from tqdm import tqdm

from typing import Callable, Iterable
from process_utils.calc import calc_acorr_order_2


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
