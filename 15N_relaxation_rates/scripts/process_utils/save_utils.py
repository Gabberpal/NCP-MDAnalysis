import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from tqdm import tqdm
from typing import Callable, Iterable, Optional
from matplotlib.backends.backend_pdf import PdfPages

from process_utils.calc import calc_acorr_order_2
from process_utils.fit import repeated_fit_auto_correlation, __multi_exp_f
from process_utils.plot import add_relpath_to_top_corner, settings_plot, get_autocorr_graph_label


def calc_and_save_acorr(path_to_vector_csv_files: Iterable[str],
                        dt_ns: float,
                        acorr_func_limit: int = -1,
                        thumbling_time: float = None,
                        acorr_func: Callable = calc_acorr_order_2,
                        out_dir: str = ".") -> None:
    """
    This function calculates the autocorrelation function for a vector array stored in CSV files and
    saves the results to an output directory.

    :param path_to_vector_csv_files: Path to a single CSV file or a list of CSV file paths containing
           the 3-dimensional vector arrays [N x 3].
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
    # Initialize index to store the length of the autocorrelation array
    index = None

    # Iterate over each CSV file containing vector data
    for path_to_vector_file in tqdm(sorted(path_to_vector_csv_files)):

        out_name = os.path.basename(path_to_vector_file).split("_")[0]
        vectors = pd.read_csv(path_to_vector_file).values
        acorr = acorr_func(vectors)[:acorr_func_limit]

        # If this is the first file, initialize the time array
        if index is None:
            index = len(acorr)
            time_ns = np.linspace(0, index * dt_ns, index, endpoint=False)

        # Reintroduce tumbling effects if thumbling_time is provided
        if (thumbling_time is not None) and (thumbling_time > 0):
            acorr *= np.exp(-time_ns / thumbling_time)

        # Create the output directory if it doesn't exist
        os.makedirs(out_dir, exist_ok=True)

        # Save the autocorrelation results to a CSV file

        pd.DataFrame({
            "time_ns": time_ns,  # Time values in nanoseconds
            "acorr": acorr  # Autocorrelation values
        }).to_csv(os.path.join(out_dir, out_name), index=False)


def fit_and_save_acorr_func(path_to_acorr_files: str,
                            bounds: list[tuple[float, float]],
                            p0: Optional[list[float]] = None,
                            lag_spacing: str = "log",
                            n_lag_points: Optional[int] = None,
                            output_directory: str = "./",
                            limit_ns: Optional[float] = None) -> None:
    """
    Fits autocorrelation functions from CSV files and saves the results.

    This function reads autocorrelation data from CSV files, fits the data using a specified
    fitting function, and saves the fitted parameters (amplitudes and relaxation times) to
    a new CSV file.

    Args:
        path_to_acorr_files (str): Path to the directory containing autocorrelation CSV files.
        bounds (List[Tuple[float, float]]): List of tuples specifying the bounds for the fitting parameters.
        p0 (Optional[List[float]]): Initial guesses for the fitting parameters. Defaults to None.
        lag_spacing (str): Spacing type for lag points. Options: "log" (logarithmic) or "linear". Defaults to "log".
        n_lag_points (Optional[int]): Number of lag points to use if lag_spacing is "log". Defaults to None.
        output_directory (str): Directory to save the output CSV files. Defaults to "./".
        limit_ns (Optional[float]): Maximum time (in nanoseconds) to consider for fitting. Defaults to None.

    Returns:
        None: The function saves the results to CSV files in the specified output directory.
    """
    # Find all CSV files in the specified directory
    path_to_ccr_csv_files = sorted(glob(os.path.join(path_to_acorr_files, "*.csv")))
    rname_list = []

    # Extract residue names from filenames
    for file in path_to_ccr_csv_files:
        aa_name = file.split('-')[2]
        rname_list.append(aa_name)

    # Iterate over bounds and fit autocorrelation functions
    for bound in bounds:
        tau_table = pd.DataFrame()
        for acorr_corr_file in tqdm(path_to_ccr_csv_files, desc=output_directory):
            df = pd.read_csv(acorr_corr_file)

            # Apply time limit if specified
            if limit_ns:
                df = df[df["time_ns"] <= limit_ns]

            time_ns, acorr = df["time_ns"].values, df["acorr"].values

            # Apply logarithmic spacing if specified
            if lag_spacing == "log":
                lag_index = np.unique(
                    np.logspace(0, int(np.log10(time_ns.size)), n_lag_points, endpoint=False).astype(int))
                acorr = np.take(acorr, lag_index)
                time_ns = np.take(time_ns, lag_index)

            # Fit the autocorrelation function
            popt = repeated_fit_auto_correlation(acorr, time_ns, bound, p0)
            name = os.path.splitext(os.path.basename(acorr_corr_file))[0]
            amplitudes = popt[::2]
            taus = popt[1::2]
            order = (len(bound[0]) + 1) // 2

            # Extract residue ID
            rid = int(name.split("-")[1])

            # Create a dictionary for the fitted parameters
            popt_dict = {
                'rId': rid,
                'rName': rname_list[rid - 2],
                'limit_ns': limit_ns
            }

            # Add amplitudes and relaxation times to the dictionary
            popt_dict.update(
                {("exp-%d-a%d" % (order, i + 1)): a for i, a in enumerate(amplitudes)}
            )
            popt_dict.update(
                {("exp-%d-tau%d" % (order, i + 1)): tau for i, tau in enumerate(taus)}
            )

            # Append the results to the DataFrame
            tau_table = pd.concat([tau_table, pd.DataFrame(popt_dict, index=[0])])

        # Sort the DataFrame by residue ID and save to CSV
        tau_table.sort_values(by=['rId'], inplace=True)
        os.makedirs(output_directory, exist_ok=True)
        tau_table.to_csv(os.path.join(output_directory, 'tau_%d_exp.csv' % order), index=False)


def plot_and_save_acorr_with_fit(path_to_fit_csv: str,
                                 path_to_acorr_csv: str,
                                 output_directory: str,
                                 ) -> None:
    """
    Reads autocorrelation data and fit parameters from CSV files, plots the autocorrelation
    with the fitted curve, and saves the plots as PDF files in the specified output directory.

    Args:
        path_to_fit_csv (str): Path to the directory containing CSV files with fit parameters.
        path_to_acorr_csv (str): Path to the directory containing CSV files with autocorrelation data.
        output_directory (str): Path to the directory where the output PDF files will be saved.

    Returns:
        None
    """
    # Create the output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)

    # Find all CSV files in the fit directory
    exp_order_acorrs_csv = glob(os.path.join(path_to_fit_csv, "*.csv"))

    # Iterate over each fit CSV file
    for tau_order_fit in sorted(exp_order_acorrs_csv):
        # Extract the base name of the file (without extension)
        name = os.path.basename(tau_order_fit).split(".")[0]

        # Create a PDF file for the current fit
        with PdfPages(os.path.join(output_directory, name + ".pdf")) as pdf:
            # Read the fit data from the CSV file
            csv_fit = os.path.join(path_to_fit_csv, name + ".csv")
            fit = pd.read_csv(csv_fit)

            # Iterate over each row in the fit data
            for (ind, (_, fit_line)) in enumerate(tqdm(fit.iterrows(), desc="plot"), 3):
                # Find the corresponding autocorrelation file
                acorr_file = glob("{}/*{:02d}*.csv".format(path_to_acorr_csv, int(fit_line["rId"])))[0]
                acorr_df = pd.read_csv(acorr_file)

                # Extract time and autocorrelation values
                time, acorr = acorr_df["time_ns"], acorr_df["acorr"]

                # Generate a label for the graph
                graph_label = get_autocorr_graph_label(fit_line)

                # Create a plot with custom settings
                fig, ax = settings_plot(graph_label)
                ax.set_title("Autocorrelation {rid} {rname}".format(
                    rid=fit_line["rId"],
                    rname=fit_line["rName"],
                ))

                # Set y-axis limits and enable grid
                ax.set_ylim(-0.1, 1.1)
                ax.grid(color="grey", alpha=0.3)

                # Extract amplitude and tau values from the fit data
                amplitude = fit_line.filter(like='-a')
                tau = fit_line.filter(like='-tau')

                # Plot the autocorrelation data and the fitted curve
                ax.plot(time, acorr)
                ax.plot(time, __multi_exp_f(time, amplitude, tau, C=0))

                # Add a vertical line at the limit_ns value
                ax.axvline(fit_line["limit_ns"], color="palegreen", ls="--")

                # Add the relative path to the top corner of the figure
                add_relpath_to_top_corner(fig)

                # Save the figure to the PDF file
                pdf.savefig(fig)
                plt.close(fig)  # Close the figure to free memory
