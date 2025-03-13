import numpy as np


def Y2m2(vectors):
    vectors = vectors / np.linalg.norm(vectors, axis=1)[:, np.newaxis]
    coef = 1.0 / 4.0 * np.sqrt(15.0 / 2.0 / np.pi)

    return coef * (vectors[:, 0] - 1j * vectors[:, 1]) ** 2


def Y2m1(vectors):
    vectors = vectors / np.linalg.norm(vectors, axis=1)[:, np.newaxis]
    coef = 1.0 / 2.0 * np.sqrt(15.0 / 2.0 / np.pi)

    return coef * (vectors[:, 0] - 1j * vectors[:, 1]) * vectors[:, 2]


def Y20(vectors):
    vectors = vectors / np.linalg.norm(vectors, axis=1)[:, np.newaxis]
    coef = 1.0 / 4.0 * np.sqrt(5.0 / np.pi)

    return coef * (2.0 * vectors[:, 2] ** 2 - vectors[:, 0] ** 2 - vectors[:, 1] ** 2)


def Y2p1(vectors):
    vectors = vectors / np.linalg.norm(vectors, axis=1)[:, np.newaxis]
    coef = -1.0 / 2.0 * np.sqrt(15.0 / 2.0 / np.pi)

    return coef * (vectors[:, 0] + 1j * vectors[:, 1]) * vectors[:, 2]


def Y2p2(vectors):
    vectors = vectors / np.linalg.norm(vectors, axis=1)[:, np.newaxis]
    coef = 1.0 / 4.0 * np.sqrt(15.0 / 2.0 / np.pi)

    return coef * (vectors[:, 0] + 1j * vectors[:, 1]) ** 2


def calc_acorr(x):
    # return real part of autocorrelation function
    f = np.fft.fft(np.pad(x, len(x), mode='constant'))
    result = np.fft.ifft(f * np.conj(f))
    result = result[:len(x)]
    result /= np.linspace(len(x), 1, len(x))
    return np.real(result)


def calc_acorr_order_2(vectors):
    vectors = vectors / np.linalg.norm(vectors, axis=1)[:, np.newaxis]
    res = [calc_acorr(vectors) for vectors in [Y2m2(vectors), Y2m1(vectors), Y20(vectors)]]
    res[0] = 2.0 * res[0]
    res[1] = 2.0 * res[1]
    return 4.0 * np.pi / 5.0 * np.sum(res, axis=0)
