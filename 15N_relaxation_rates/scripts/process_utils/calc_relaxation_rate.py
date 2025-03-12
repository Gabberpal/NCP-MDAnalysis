import os
import pandas as pd
import numpy as np
from math import pi


def _calc_R1(amplitude, taus_s, nmr_freq):
    """
    Calculates the longitudinal relaxation rate (R1) for NMR spectroscopy based on spectral density functions.

    Args:
        amplitude (list or np.array): Amplitudes of the spectral density components.
        taus_s (list or np.array): Correlation times (in seconds) for the spectral density components.
        nmr_freq (float): NMR frequency in Hz.

    Returns:
        float: The calculated R1 relaxation rate.
    """

    def J(w):
        """
        Calculates the spectral density function J(w) for a given frequency w.

        Args:
            w (float): Frequency in rad/s.

        Returns:
            float: The value of the spectral density function at frequency w.
        """
        return sum(
            (2 / 5 * a * tau / (1 + np.square(tau * w))) for a, tau in zip(amplitude, taus_s)
        )

    # Physical constants
    u0 = 4 * pi * 1e-7  # Vacuum permeability (H/m)
    h = 6.626069e-34  # Planck's constant (J·s)
    gH = 267.522e6  # Gyromagnetic ratio of hydrogen (rad/(s·T))
    gX = -27.126e6  # Gyromagnetic ratio of the heteronucleus (rad/(s·T))
    CSA = -172.0  # Chemical shift anisotropy (ppm)
    rHX = 1.02e-10  # Distance between hydrogen and heteronucleus (m)

    # Calculate the dipolar coupling constant (d^2)
    d2 = ((u0 / 4 / pi) * (h / 2 / pi) * (gH * gX) / (rHX ** 3)) ** 2

    # Calculate the ratio of gyromagnetic ratios (gamma_H / gamma_X)
    gHX = gH / gX

    # Calculate angular frequencies for hydrogen (wH) and heteronucleus (wX)
    wH = 2 * pi * nmr_freq  # Angular frequency of hydrogen (rad/s)
    wX = wH / gHX  # Angular frequency of heteronucleus (rad/s)

    # Calculate the chemical shift anisotropy term (c^2)
    c2 = (1.0 / 3.0) * ((CSA * 1e-6 * wX) ** 2)

    # Calculate the R1 relaxation rate using spectral density functions
    R1 = 0.25 * d2 * (3 * J(wX) + J(wH - wX) + 6 * J(wH + wX)) + c2 * J(wX)

    return R1


def _calc_R2(amplitude, taus_s, nmr_freq):
    """
    Calculates the transverse relaxation rate (R2) for NMR spectroscopy based on spectral density functions.

    Args:
        amplitude (list, np.array, int, or float): Amplitudes of the spectral density components.
        taus_s (list, np.array, int, or float): Correlation times (in seconds) for the spectral density components.
        nmr_freq (float): NMR frequency in Hz.

    Returns:
        float: The calculated R2 relaxation rate.
    """

    def J(w):
        """
        Calculates the spectral density function J(w) for a given frequency w.

        Args:
            w (float): Frequency in rad/s.

        Returns:
            float: The value of the spectral density function at frequency w.
        """
        return sum(
            (2 / 5 * a * tau / (1 + np.square(tau * w))) for a, tau in zip(amplitude, taus_s)
        )

    # Physical constants
    u0 = 4 * pi * 1e-7  # Vacuum permeability (H/m)
    h = 6.626069e-34  # Planck's constant (J·s)
    gH = 267.522e6  # Gyromagnetic ratio of hydrogen (rad/(s·T))
    gX = -27.126e6  # Gyromagnetic ratio of the heteronucleus (rad/(s·T))
    CSA = -172.0  # Chemical shift anisotropy (ppm)
    rHX = 1.02e-10  # Distance between hydrogen and heteronucleus (m)

    # Calculate the dipolar coupling constant (d^2)
    d2 = ((u0 / 4 / pi) * (h / 2 / pi) * (gH * gX) / (rHX ** 3)) ** 2

    # Calculate the ratio of gyromagnetic ratios (gamma_H / gamma_X)
    gHX = gH / gX

    # Calculate angular frequencies for hydrogen (wH) and heteronucleus (wX)
    wH = 2 * pi * nmr_freq  # Angular frequency of hydrogen (rad/s)
    wX = wH / gHX  # Angular frequency of heteronucleus (rad/s)

    # Calculate the chemical shift anisotropy term (c^2)
    c2 = (1.0 / 3.0) * ((CSA * 1e-6 * wX) ** 2)

    # Ensure amplitude and taus_s are iterable (convert to lists if they are single values)
    if isinstance(amplitude, (int, float)):
        amplitude = [amplitude]
    if isinstance(taus_s, (int, float)):
        taus_s = [taus_s]

    # Calculate the R2 relaxation rate using spectral density functions
    R2 = 0.125 * d2 * (4 * J(0) + 3 * J(wX) + J(wH - wX) + 6 * J(wH) + 6 * J(wH + wX)) + (1.0 / 6.0) * c2 * (
                4 * J(0) + 3 * J(wX))

    return R2


def calc_relaxition(path_to_fit_csv, nmr_freq, func):
    """
    Calculates relaxation rates (R1 or R2) for multiple entries in a CSV file containing fit parameters.

    Args:
        path_to_fit_csv (str): Path to the CSV file containing fit parameters (amplitudes and taus).
        nmr_freq (float): NMR frequency in Hz.
        func (callable): Function to calculate the relaxation rate (e.g., `_calc_R1` or `_calc_R2`).

    Returns:
        pd.DataFrame: A DataFrame containing the calculated relaxation rates for each entry in the CSV file.
                      Columns include 'rName', 'rId', and 'relaxation_rate'.
    """
    # Read the CSV file containing fit parameters
    df = pd.read_csv(path_to_fit_csv)

    # Initialize an empty DataFrame to store the results
    rate_table = pd.DataFrame()

    # Iterate over each row in the DataFrame
    for ind, fit_line in df.iterrows():
        # Extract amplitudes and taus from the fit line
        amplitude = fit_line.filter(like='-a')  # Amplitudes (columns containing '-a')
        taus = fit_line.filter(like='-tau')  # Correlation times (columns containing '-tau')
        taus_s = taus * 1e-9  # Convert taus from nanoseconds to seconds

        # Calculate the relaxation rate using the provided function
        rate = func(amplitude, taus_s, nmr_freq)

        # Create a dictionary with the results for the current entry
        D = {
            'rName': df.rName[ind],  # Name of the residue
            'rId': df.rId[ind],  # ID of the residue
            "relaxation_rate": rate  # Calculated relaxation rate
        }

        # Append the results to the rate_table DataFrame
        temp = pd.DataFrame(D, index=[0])
        rate_table = pd.concat([rate_table, temp])

    return rate_table


def get_relaxition_rate(path_to_fit, nmr_freq, output_directory="./", rate="R1"):
    """
    Calculates relaxation rates (R1 or R2) for multiple fit files and saves the results to a CSV file.

    Args:
        path_to_fit (str): Path to the directory containing fit files (e.g., `tau_*_exp.csv`).
        nmr_freq (float): NMR frequency in Hz.
        output_directory (str, optional): Directory to save the output CSV file. Defaults to "./".
        rate (str, optional): Type of relaxation rate to calculate. Must be "R1" or "R2". Defaults to "R1".

    Returns:
        pd.DataFrame: A DataFrame containing the calculated relaxation rates for all fit files.
    """
    from glob import glob

    # Dictionary mapping rate types to their corresponding calculation functions
    func_dict = {"R1": _calc_R1, "R2": _calc_R2}

    # Initialize an empty DataFrame to store the results
    df = pd.DataFrame()

    # Define the pattern for fit files
    fits = ["tau_*_exp.csv"]

    # Iterate over each fit file pattern
    for fit in fits:
        # Construct the full path to the fit files
        path_to_fit_csv = os.path.join(path_to_fit, fit)

        # Find the most recent fit file matching the pattern
        latest_fit_file = sorted(glob(path_to_fit_csv))[-1]

        # Calculate relaxation rates using the specified function
        relaxation_rate = calc_relaxition(latest_fit_file, nmr_freq, func_dict[rate])

        # Merge the results into the main DataFrame
        if df.empty:
            df = relaxation_rate
        else:
            df = pd.merge(df, relaxation_rate, left_index=False, right_index=False)

    # Create the output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)

    # Save the results to a CSV file
    output_file = os.path.join(output_directory, "{rate}.csv".format(rate=rate))
    df.to_csv(output_file, index=False)

    return df


if __name__ == '__main__':
    nmr_freq = 850e6
    path_to_fit = "data/fit/"
    for rate in ["R1", "R2"]:
        get_relaxition_rate(path_to_fit=path_to_fit, nmr_freq=nmr_freq, rate=rate)
