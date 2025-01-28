import numpy as np
from scipy.optimize import curve_fit
from typing import Tuple, List, Union, Iterable


def __multi_exp_f(x: Union[float, int],
                  A: List[Union[float, int]],
                  TAU: List[Union[float, int]],
                  C: Union[float, int]) -> float:
    """
    :param x: argument of some exponential functions composition
    :param A: array of amplitudes
    :param TAU: array of time constants
    :param C: free element
    :return: sum exponential functions composition
    """
    return sum(
        (a * np.exp(-x / tau)) for a, tau in zip(A, TAU)
    ) + C


def multi_exp_fixed_amplitude_1(x: Union[float, int], *args):
    """
    :param x: argument of some exponential functions composition
    :param args: array of amplitudes and time constants
    :return: callable __multi_exp_f
    """
    TAU = args[0::2]

    if len(args) % 2 == 1:
        C = 0
        A = args[1:-1:2]
    else:
        C = args[-1]
        A = args[1:-1:2]
    A0 = 1 - sum(A) - C
    return __multi_exp_f(x, [A0] + list(A), TAU, C)


def fit_auto_correlation(time: List[float],
                         acorr: List[float],
                         bounds: List[List[List[Union[float, int]]]],
                         p0=None
                         ) \
        -> Tuple[int, Union[np.ndarray, Iterable, int, float]]:
    """
    Fit input data with :math:`\\sum_n A_n \\exp(-t/\\tau_n) + const`

    :param time: time data series
    :param acorr: auto-correlation data series
    :param bounds: curve parameters bounds
    :return: Fit curve parameters
    """

    if p0 is None:
        p0 = np.mean(bounds, axis=0)[1:]

    args, pcov = curve_fit(multi_exp_fixed_amplitude_1,
                           time,
                           acorr,
                           p0=p0,
                           bounds=np.array(bounds)[:, 1:],
                           max_nfev=10000,
                           ftol=1e-6)

    if len(args) % 2 == 1:
        C = 0
        A = args[1:-1:2]
    else:
        C = args[-1]
        A = args[1:-1:2]
    A0 = 1 - sum(A) - C

    return [A0] + list(args)


def repeated_fit_auto_correlation(corr: List[float],
                                  time: List[float],
                                  bounds: List[List[Union[float, int]]],
                                  p0,
                                  repeats=3,
                                  fit_func=multi_exp_fixed_amplitude_1,
                                  ) \
        -> Tuple[int, Union[np.ndarray, Iterable, int, float]]:
    """
    Fit input data with :math:`\\sum_n A_n \\exp(-t/\\tau_n) + const`

    :param time: time data series
    :param corr: correlation data series
    :param bounds: curve parameters bounds
    :param repeats: number of fit repetition
    :repeated_fit_auto_correlation: function for repeated fit data with initial guess subjected to a random variation
    :return: Fit curve parameters
    """

    def randminmax(min, max):
        return np.random.rand(1) * (max - min) + min

    R_square = []
    popt_all = []

    if p0 is None:
        p0 = np.mean(bounds, axis=0)[1:]

    for _ in range(repeats):

        try:
            popt = fit_auto_correlation(time, corr, bounds, p0)
            R_square.append(sum((np.array(corr) - np.array(fit_func(time, *popt))) ** 2))
            popt_all.append(popt)
            p0[1::2][0] = randminmax(p0[1::2][0] / 10, p0[1::2][0])
            p0[1::2][1:] = randminmax(p0[1::2][:-1], p0[1::2][1:])

        except RuntimeError:
            print("Fit error n={}".format(len(bounds[0]) // 2))

    min_ind_r_square = np.argmin(R_square)

    return popt_all[min_ind_r_square]
