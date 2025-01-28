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
    d2 = ((u0 / 4 / pi) * (h / 2 / pi) * (gH * gX) / (rHX**3)) ** 2

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
    d2 = ((u0 / 4 / pi) * (h / 2 / pi) * (gH * gX) / (rHX**3)) ** 2

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
    R2 = 0.125 * d2 * (4 * J(0) + 3 * J(wX) + J(wH - wX) + 6 * J(wH) + 6 * J(wH + wX)) + \
         (1.0 / 6.0) * c2 * (4 * J(0) + 3 * J(wX))

    return R2


def calc_relaxition(path_to_fit_csv, nmr_freq, func):

    df = pd.read_csv(path_to_fit_csv)
    rate_table = pd.DataFrame()
    for ind,fit_line in df.iterrows():
        amplitude = fit_line.filter(like='-a') 
        taus = fit_line.filter(like='-tau')
        taus_s = taus * 1e-9
        rate = func(amplitude, taus_s, nmr_freq)
        D = {'rName': df.rName[ind], 'rId': df.rId[ind], "relaxation_rate": rate}
        temp = pd.DataFrame(D, index=[0])
        rate_table = pd.concat([rate_table, temp])
    return rate_table


def get_relaxition_rate(path_to_fit, nmr_freq, output_directory="./", rate="R1"):
    from glob import glob
    func_dict = {"R1": _calc_R1, "R2": _calc_R2}
    df = pd.DataFrame()
    fits = ["tau_*_exp.csv"]
    for fit in fits:
        path_to_fit_csv = os.path.join(path_to_fit, fit)
        relaxation_rate = calc_relaxition(sorted(glob(path_to_fit_csv))[-1], nmr_freq, func_dict[rate])
        if df.empty:
            df = relaxation_rate
        else:
            df = pd.merge(df, relaxation_rate, left_index=False, right_index=False)
    os.makedirs(output_directory, exist_ok=True)
    df.to_csv(os.path.join(output_directory, "{rate}.csv".format(rate=rate)), index=False)
    return df



if __name__ == '__main__':
    nmr_freq = 850e6
    path_to_fit = "data/fit/"
    for rate in ["R1", "R2"]:
        get_relaxition_rate(path_to_fit=path_to_fit, nmr_freq=nmr_freq, rate=rate)
