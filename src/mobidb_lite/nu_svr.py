# Giulio Tesei, Sören von Bülow, Kresten Lindorff-Larsen 2024

from joblib import dump, load
import numpy as np
import pandas as pd
import random
import numba as nb


@nb.jit(nopython=True)
def calc_fcr_and_ncpr(qs):
    fcr = np.mean(np.abs(qs))
    ncpr = np.mean(qs)
    return fcr, ncpr


@nb.jit(nopython=True)
def calc_sigma(qs):
    fcr, ncpr = calc_fcr_and_ncpr(qs)
    if fcr == 0:
        return 0.
    else:
        return ncpr ** 2 / fcr


@nb.jit(nopython=True)
def calc_delta_form(qs, window=5):
    sig_m = calc_sigma(qs)

    nw = len(qs) - window + 1
    sigs = np.zeros((nw))

    for idx in range(0, nw):
        q_window = qs[idx:idx + window]
        sigs[idx] = calc_sigma(q_window)
    delta = np.sum((sigs - sig_m) ** 2) / nw
    return delta


@nb.jit(nopython=True)
def calc_delta(qs):
    d5 = calc_delta_form(qs, window=5)
    d6 = calc_delta_form(qs, window=6)
    return (d5 + d6) / 2.


@nb.jit(nopython=True)
def get_qs_fast(seq):
    qs = np.zeros(len(seq))

    # loop through sequence
    for idx in range(len(seq)):
        if seq[idx] in ['R', 'K']:
            qs[idx] = 1.
        elif seq[idx] in ['E', 'D']:
            qs[idx] = -1.
        else:
            qs[idx] = 0.
    return qs


@nb.jit(nopython=True)
def check_dmax(seq, dmax, seqmax):
    qs = get_qs_fast(seq)
    d = calc_delta(qs)
    if d > dmax:
        return seq, d
    else:
        return seqmax, dmax


@nb.jit(nopython=True)
def calc_case0(seqpos, seqneg, seqneu):
    seqmax = ''
    dmax = 0.
    N = len(seqpos) + len(seqneg) + len(seqneu)
    if len(seqpos) == 0:
        seqcharge = seqneg[:]
    elif len(seqneg) == 0:
        seqcharge = seqpos[:]
    if len(seqneu) > len(seqcharge):
        for pos in range(0, N - len(seqcharge) + 1):
            seqout = ''
            seqout += seqneu[:pos]
            seqout += seqcharge
            seqout += seqneu[pos:]
            seqmax, dmax = check_dmax(seqout, dmax, seqmax)
    else:
        for pos in range(0, N - len(seqneu) + 1):
            seqout = ''
            seqout += seqcharge[:pos]
            seqout += seqneu
            seqout += seqcharge[pos:]
            seqmax, dmax = check_dmax(seqout, dmax, seqmax)
    return seqmax


@nb.jit(nopython=True)
def calc_case1(seqpos, seqneg, seqneu):
    seqmax = ''
    dmax = 0.
    N = len(seqpos) + len(seqneg) + len(seqneu)
    if len(seqpos) > len(seqneg):
        for pos in range(0, N - len(seqneg) + 1):
            seqout = ''
            seqout += seqpos[:pos]
            seqout += seqneg
            seqout += seqpos[pos:]
            seqmax, dmax = check_dmax(seqout, dmax, seqmax)
    else:
        for neg in range(0, N - len(seqpos) + 1):
            seqout = ''
            seqout += seqneg[:neg]
            seqout += seqpos
            seqout += seqneg[neg:]
            seqmax, dmax = check_dmax(seqout, dmax, seqmax)
    return seqmax


@nb.jit(nopython=True)
def calc_case2(seqpos, seqneg, seqneu):
    seqmax = ''
    dmax = 0.
    for startNeuts in range(0, 7):
        for endNeuts in range(0, 7):
            startBlock = seqneu[:startNeuts]
            endBlock = seqneu[startNeuts:startNeuts + endNeuts]
            midBlock = seqneu[startNeuts + endNeuts:]

            seqout = ''
            seqout += startBlock
            seqout += seqpos
            seqout += midBlock
            seqout += seqneg
            seqout += endBlock
            seqmax, dmax = check_dmax(seqout, dmax, seqmax)
    return seqmax


@nb.jit(nopython=True)
def calc_case3(seqpos, seqneg, seqneu):
    seqmax = ''
    dmax = 0.
    for midNeuts in range(0, len(seqneu) + 1):
        midBlock = seqneu[:midNeuts]
        for startNeuts in range(0, len(seqneu) - midNeuts + 1):
            startBlock = seqneu[midNeuts:midNeuts + startNeuts]
            seqout = ''
            seqout += startBlock
            seqout += seqpos
            seqout += midBlock
            seqout += seqneg
            seqout += seqneu[midNeuts + startNeuts:]
            seqmax, dmax = check_dmax(seqout, dmax, seqmax)
    return seqmax


def shuffle_str(seq):
    l = list(seq)
    random.shuffle(l)
    seq_shuffled = "".join(l)
    return (seq_shuffled)


def split_seq(seq):
    """ split sequence in positive, negative, neutral residues """
    seqpos = []
    seqneg = []
    seqneu = []
    for s in seq:
        if s in ['K', 'R']:
            seqpos.append(s)
        elif s in ['D', 'E']:
            seqneg.append(s)
        else:
            seqneu.append(s)
    seqpos = shuffle_str(seqpos)
    seqneg = shuffle_str(seqneg)
    seqneu = shuffle_str(seqneu)
    return seqpos, seqneg, seqneu


def construct_deltamax(seq):
    seqpos, seqneg, seqneu = split_seq(seq)

    if (len(seqpos) == 0) or (len(seqneg) == 0):
        seqmax = calc_case0(seqpos, seqneg, seqneu)
    elif len(seqneu) == 0:
        seqmax = calc_case1(seqpos, seqneg, seqneu)
    elif len(seqneu) >= 18:
        seqmax = calc_case2(seqpos, seqneg, seqneu)
    else:
        seqmax = calc_case3(seqpos, seqneg, seqneu)
    return seqmax


def calc_fcr_and_kappa(seq):
    """ kappa parameter, Das & Pappu, PNAS 2013
    Here set to 0, instead of -1 as in localCIDER, for sequences
    with no charges or only pos/neg residues """
    qs = get_qs_fast(seq)
    if np.all(qs == 0):
        return 0, 0
    else:
        fcr, _ = calc_fcr_and_ncpr(qs)
        seqpos, seqneg, seqneu = split_seq(seq)
        if (len(seqneu) == 0):
            if (len(seqneg) == 0) or (len(seqpos) == 0):
                return fcr, 0

    delta = calc_delta(qs)

    seq_max = construct_deltamax(seq)
    qs_max = get_qs_fast(seq_max)
    delta_max = calc_delta(qs_max)

    kappa = delta / delta_max
    return fcr, kappa


@nb.jit(nopython=True)
def calc_SCD(seq):
    """ Sequence charge decoration, eq. 14 in Sawle & Ghosh, JCP 2015 """
    qs = get_qs_fast(seq)
    N = len(seq)
    scd = 0.
    for idx in range(1, N):
        for jdx in range(0, idx):
            s = qs[idx] * qs[jdx] * (idx - jdx) ** 0.5
            scd = scd + s
    scd = scd / N
    return scd


@nb.jit(nopython=True)
def calc_SHD(seq, lambda_map, beta=-1.):
    """ Sequence hydropathy decoration, eq. 4 in Zheng et al., JPC Letters 2020"""
    N = len(seq)
    shd = 0.
    for idx in range(0, N - 1):
        for jdx in range(idx + 1, N):
            s = lambda_map[idx, jdx] * (jdx - idx) ** beta
            shd = shd + s
    shd = shd / N
    return shd


def calc_seq_properties(seq):
    lambdas = {"R": 0.7307624768, "D": 0.0416040481, "N": 0.425585901, "E": 0.0006935461,
               "K": 0.1790211739, "H": 0.4663667291, "Q": 0.3934318551, "S": 0.4625416812,
               "C": 0.5615435099, "G": 0.7058843734, "T": 0.3713162976, "A": 0.2743297969,
               "M": 0.5308481134, "Y": 0.9774611449, "V": 0.2083769608, "W": 0.989376474,
               "L": 0.6440005008, "I": 0.5423623611, "P": 0.3593126576, "F": 0.8672358982}
    lambdas = np.asarray(list(map(lambda x: lambdas[x], seq)))
    lambda_map = np.add.outer(lambdas, lambdas)
    scd = calc_SCD(seq)  # Slow
    shd = calc_SHD(seq, lambda_map)
    fcr, kappa = calc_fcr_and_kappa(seq)
    mean_lambda = np.mean(lambdas)
    return scd, shd, kappa, fcr, mean_lambda


if __name__ == "__main__":

    # load scikit-learn SVR model
    # model_nu = load('data/svr_model_nu.joblib')
    # dump(model_nu, open('data/svr_model_nu_2.joblib', 'wb'))

    model_nu = load('data/svr_model_nu_2.joblib')

    # sequence to analyse
    seq = "ISSENSNPEQDLKLTSEEESQRLKVSENSQPEKMSQEPEINKDCDREVEEEIKKHGSNPVGLPENLTNGASAGNGDDGLIPQRKSRKPENQQFPDTENEEYHSDEQNDTQKQLSEEQNTGISQDEILTNKQ"

    if len(seq) < 30 or len(seq) > 1500:
        print('Sequence length outside the range of the training set')

    else:
        # calculate sequence features
        scd, shd, kappa, fcr, mean_lambda = calc_seq_properties(seq)

        # calculate scaling exponent
        nu_svr = model_nu.predict([[scd, shd, kappa, fcr, mean_lambda]])[0]

        # print results
        results = pd.Series(data=[scd, shd, kappa, fcr, mean_lambda, nu_svr],
                            index=['scd', 'shd', 'kappa', 'fcr', 'mean_lambda', 'nu_svr'])

        print(results)
