import os
import pandas as pd
import numpy as np
from math import pi


def _calc_R1(amplitude, taus_s, nmr_freq):

    def J(w):
        return sum(
            (2 / 5 * a * tau/(1+np.square(tau*w))) for a, tau in zip(amplitude, taus_s))

    u0 = 4 * pi * 1e-7
    h = 6.626069e-34
    gH = 267.522e6
    gX = -27.126e6
    CSA=-172.0  # unit: ppm
    rHX=1.02e-10

    d2 = ((u0/4/pi)*(h/2/pi)*(gH*gX)/(rHX**3))**2
    gHX = gH / gX  # ratio gamma_H to gamma_X
    wH = 2 * pi * nmr_freq
    wX = wH/gHX
    c2 = (1./3.)*((CSA*1e-6*wX)**2)

    R1 = 0.25*d2*(3*J(wX)+J(wH-wX)+6*J(wH+wX)) + c2*J(wX)
    return R1

def _calc_R2(amplitude, taus_s, nmr_freq):
    def J(w):
        return sum(
            (2 / 5 * a * tau / (1 + np.square(tau * w))) for a, tau in zip(amplitude, taus_s))

    # Calculate dipolar coupling constant d2 -------------------------------
    u0 = 4 * pi * 1e-7
    h = 6.626069e-34
    gH = 267.522e6
    gX = -27.126e6
    CSA = -172.0  # unit: ppm
    rHX = 1.02e-10

    d2 = ((u0 / 4 / pi) * (h / 2 / pi) * (gH * gX) / (rHX ** 3)) ** 2
    gHX = gH / gX  # ratio gamma_H to gamma_X
    wH = 2 * pi * nmr_freq
    wX = wH / gHX
    c2 = (1. / 3.) * ((CSA * 1e-6 * wX) ** 2)

    if isinstance(amplitude, int) or isinstance(amplitude, float):
        amplitude = [amplitude]
    if isinstance(taus_s, int) or isinstance(taus_s, float):
        taus_s = [taus_s]

    R2 = 0.125 * d2 * (4 * J(0) + 3 * J(wX) + J(wH - wX) + 6 * J(wH) + 6 * J(wH + wX)) + \
         (1. / 6.) * c2 * (4 * J(0) + 3 * J(wX))
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
