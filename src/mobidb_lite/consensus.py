import math
import os
import re
import sys
from tempfile import mkstemp
from typing import Union
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import reduce

from mobidb_lite import disembl, espritz, globplot, iupred, seg, anchor

try:
    from mobidb_lite import nu_svr
except ImportError:
    print("Nu svr module not loaded")
    pass


_THRESHOLDS = {
    "mobidblite": 0.625,
    "disembl-hl": 0.086,
    "disembl-rem465": 0.5,
    "espritz-d": 0.5072,
    "espritz-n": 0.3089,
    "espritz-x": 0.1434,
    "globplot": 0,
    "iupred-l": 0.5,
    "iupred-s": 0.5,
    "seg": 0.5,
    "anchor": 0.5
}

_MOBIDB_NAMES = {
    "majority": "prediction-disorder-th_50",
    "mobidblite": "prediction-disorder-mobidb_lite",
    "disembl-hl": "prediction-disorder-disHL",
    "disembl-rem465": "prediction-disorder-dis465",
    "espritz-d": "prediction-disorder-espD",
    "espritz-n": "prediction-disorder-espN",
    "espritz-x": "prediction-disorder-espX",
    "globplot": "prediction-disorder-glo",
    "iupred-l": "prediction-disorder-iupl",
    "iupred-s": "prediction-disorder-iups",
    "seg": "prediction-low_complexity-seg",
    "anchor": "prediction-lip-anchor",
    "Polyampholyte": "prediction-polyampholyte-mobidb_lite_sub",
    "Positive Polyelectrolyte": "prediction-positive_polyelectrolyte-mobidb_lite_sub",
    "Negative Polyelectrolyte": "prediction-negative_polyelectrolyte-mobidb_lite_sub",
    "Cysteine-rich": "prediction-cysteine_rich-mobidb_lite_sub",
    "Proline-rich": "prediction-proline_rich-mobidb_lite_sub",
    "Glycine-rich": "prediction-glycine_rich-mobidb_lite_sub",
    "Low complexity": "prediction-low_complexity-mobidb_lite_sub",
    "Polar": "prediction-polar-mobidb_lite_sub",
    "Extended": "prediction-extended-mobidb_lite_sub",
    "Compact": "prediction-compact-mobidb_lite_sub"
}

_CONSENSUS = frozenset(
    ["disembl-hl", "disembl-rem465", "espritz-d", "espritz-n", "espritz-x", "globplot", "iupred-l", "iupred-s"])

_POSITIVE_FLAG = "1"
_NEGATIVE_FLAG = "0"
_POSITIVE_RESIDUES = {"H", "R", "K"}
_NEGATIVE_RESIDUES = {"D", "E"}
_FEATURES = [
    "Polyampholyte",
    "Positive Polyelectrolyte",
    "Negative Polyelectrolyte",
    "Cysteine-rich",
    "Proline-rich",
    "Glycine-rich",
    "Low complexity",
    "Polar",
    "Extended",
    "Compact"
]

_VALID_AA = re.compile('[^RDNEKHQSCGTAMYVWLIPF]')

def parse_fasta(file):
    seq_id = sequence = ""

    with file as fh:
        for line in map(str.rstrip, fh):
            if line[0] == ">":
                if seq_id and sequence:
                    yield seq_id, sequence.upper()
                seq_id = line[1:].split()[0]
                sequence = ""
            elif line:
                sequence += line

    if seq_id and sequence:
        yield seq_id, sequence.upper()


def run_predictors(sequence: str, bindir: str, **kwargs) -> dict:
    tempdir = kwargs.get("tempdir")
    run_extra = kwargs.get("extra", False)
    round_score = kwargs.get("round", False)

    fd, disbin = mkstemp(dir=tempdir)
    with open(fd, "wt") as fh:
        fh.write(f"1\n{len(sequence)}\n{sequence}")

    hot_loop, remark_465 = disembl.run(os.path.join(bindir, "DisEMBL"),
                                       os.path.join(bindir, "TISEAN"),
                                       sequence)

    results = {"disembl-hl": hot_loop, "disembl-rem465": remark_465,
               "espritz-d": espritz.run_espritz_d(os.path.join(bindir, "ESpritz"), disbin),
               "espritz-n": espritz.run_espritz_n(os.path.join(bindir, "ESpritz"), disbin),
               "espritz-x": espritz.run_espritz_x(os.path.join(bindir, "ESpritz"), disbin),
               "globplot": globplot.run(os.path.join(bindir, "TISEAN"), sequence),
               "iupred-l": iupred.run_long(os.path.join(bindir, "IUPred"), sequence),
               "iupred-s": iupred.run_short(os.path.join(bindir, "IUPred"), sequence)}
    # "fess": fess.run_fess(os.path.join(bindir, "FeSS"), disbin)}

    os.unlink(disbin)


    fd, fasta = mkstemp(dir=tempdir)
    with open(fd, "wt") as fh:
        fh.write(f">1\n{sequence}\n")

    results["seg"] = seg.run(os.path.join(bindir, "SEG"), fasta)

    if run_extra:
            results['anchor'] = anchor.run_anchor(os.path.join(bindir, "ANCHOR"), fasta)

    os.unlink(fasta)

    states = {}
    for k in results:
        states[k] = [(round(s, 3) if round_score else s) >= _THRESHOLDS[k] for s in results[k]]

    return results, states


def run(file: str, bindir: str, datadir: str, threads: int, **kwargs):

    # Import model nu data fro compact/extended calculation
    if kwargs.get("extra", False):
        model_nu = nu_svr.load(os.path.join(datadir, "svr_model_nu_2.joblib"))
        kwargs["model_nu"] = model_nu

    if threads > 1:
        with ThreadPoolExecutor(max_workers=threads) as executor:
            fs = {}
            for seq_id, sequence in parse_fasta(file):
                f = executor.submit(predict, seq_id, sequence, bindir, **kwargs)
                fs[f] = seq_id

                if len(fs) == 1000:
                    for f in as_completed(fs):
                        seq_id = fs[f]
                        regions, scores = f.result()
                        yield seq_id, regions, scores

                    fs.clear()

            for f in as_completed(fs):
                seq_id = fs[f]
                regions, scores = f.result()
                yield seq_id, regions, scores
    else:
        for seq_id, sequence in parse_fasta(file):
            regions, scores = predict(seq_id, sequence, bindir, **kwargs)
            yield seq_id, regions, scores


def predict(sequence_id: str, sequence: str, bindir: str, **kwargs):
    force_consensus = kwargs.get("force", False)
    round_score = kwargs.get("round", False)
    run_extra = kwargs.get("run_extra", False)
    tempdir = kwargs.get("tempdir")
    threshold = kwargs.get("threshold", _THRESHOLDS["mobidblite"])
    model_nu = kwargs.get("model_nu", None)

    seq_length = len(sequence)
    scores, scores_states = run_predictors(sequence, bindir, extra=run_extra, tempdir=tempdir)

    # Sub-regions
    if len(sequence) == len(scores_states["seg"]):
        features = get_region_features(sequence, scores_states["seg"])
    else:
        features = None

    # Predictors agreement (subset of methods)
    agreement = [0] * seq_length
    num_indicators = 0

    for pred_name, pred_state in scores_states.items():
        if pred_name in _CONSENSUS:
            if pred_state is None or len(pred_state) != seq_length:
                sys.stderr.write(f"{sequence_id}: {pred_name} excluded\n")
                continue

            num_indicators += 1
            for i, state in enumerate(pred_state):
                agreement[i] += state

    if num_indicators == 0:
        return None, None
    elif num_indicators < len(scores) and not force_consensus:
        return None, None

    # Consensus states
    states = ""
    states_majority = ""
    scores["mobidblite"] = []
    for s in agreement:
        if round_score:
            score = round(s / num_indicators, 3)
        else:
            score = s / num_indicators
        scores["mobidblite"].append(score)

        if score >= threshold:
            states += _POSITIVE_FLAG
        else:
            states += _NEGATIVE_FLAG

        if score >= 0.5:
            states_majority += _POSITIVE_FLAG
        else:
            states_majority += _NEGATIVE_FLAG


    # Majority consensus
    scores_states["majority"] = states_majority

    # Post-processing
    states = dilate(states, max_length=3)
    states = erode(states, max_length=3)
    states = merge_long_disordered_regions(states)

    # Convert states lists into strings
    scores_states = {pred_name: "".join([str(int(s)) for s in state]) for pred_name, state in scores_states.items()}
    scores_states["mobidblite"] = states

    results = {}

    for pred_name, pred_state in scores_states.items():
        if pred_state is None or len(pred_state) != seq_length:
            sys.stderr.write(f"{sequence_id}: {pred_name} excluded\n")
            continue

        # Regions
        if pred_name == "mobidblite":
            regions = get_regions(pred_state, min_length=20)
        else:
            regions = get_regions(pred_state, min_length=1)

        for start, end, _ in sorted(regions):
            results.setdefault(pred_name, []).append((start + 1, end + 1))

            # Sub-regions
            if features and pred_name == "mobidblite":
                region = features[start:end + 1]

                for i, j, state in get_regions(region, min_length=10):
                    # state = _FEATURES[int(x)-1]
                    results.setdefault(state, []).append((start + 1 + i, start + 1 + j))

            # Nu parameters
            seq = _VALID_AA.sub('', sequence[start:end + 1])  # Remove non-standard AA from the sequence fragment
            if pred_name == "mobidblite" and model_nu and 30 <= len(seq) <= 1500:

                # calculate sequence features
                scd, shd, kappa, fcr, mean_lambda = nu_svr.calc_seq_properties(seq)
                # calculate scaling exponent
                nu = model_nu.predict([[scd, shd, kappa, fcr, mean_lambda]])[0]

                #  "compact" se nu <= 0.475 e "extended" se nu > 0.55
                if nu <= 0.475:
                    results.setdefault("Compact", []).append((start + 1, end + 1, nu))
                elif nu > 0.55:
                    results.setdefault("Extended", []).append((start + 1, end + 1, nu))

    return results, scores


def dilate(states: str, max_length: int) -> str:
    states = "{0}{1}{0}".format(_POSITIVE_FLAG * max_length, states)

    for level in range(1, max_length + 1):
        old = "{0}{1}{0}".format(_POSITIVE_FLAG * level,
                                 _NEGATIVE_FLAG * level)
        new = "{0}{0}{0}".format(_POSITIVE_FLAG * level)

        for _ in range(level + 1):
            states = states.replace(old, new)

    return states[max_length:-max_length]


def erode(states: str, max_length: int) -> str:
    states = "{0}{1}{0}".format(_NEGATIVE_FLAG * max_length, states)

    for level in range(1, max_length + 1):
        old = "{0}{1}{0}".format(_NEGATIVE_FLAG * level,
                                 _POSITIVE_FLAG * level)
        new = "{0}{0}{0}".format(_NEGATIVE_FLAG * level)

        for _ in range(level + 1):
            states = states.replace(old, new)

    return states[max_length:-max_length]


def merge_long_disordered_regions(states: str) -> str:
    while True:
        new_states = re.sub(
            pattern=r"{p}{{21,}}{n}{{1,10}}{p}{{21,}}".format(
                p=_POSITIVE_FLAG,
                n=_NEGATIVE_FLAG
            ),
            repl=_repl_struct_by_disord,
            string=states)

        if new_states == states:
            return new_states

        states = new_states


def _repl_struct_by_disord(match: re.match):
    return _POSITIVE_FLAG * len(match.group(0))


def get_regions(states: Union[list, str], min_length: int) -> list:
    # Return 0-indexed regions
    regions = []
    start = None
    current_flag = None
    for i, flag in enumerate(states):
        if flag != current_flag:
            if start is not None and current_flag != _NEGATIVE_FLAG:
                end = i - 1
                length = end - start + 1
                if length >= min_length:
                    regions.append((start, end, current_flag))

            start = i
            current_flag = flag

    if start is not None and current_flag != _NEGATIVE_FLAG:
        end = len(states) - 1
        length = end - start + 1
        if length >= min_length:
            regions.append((start, end, current_flag))

    return regions


def get_region_features(sequence: str, seg_states: list) -> list:

    all_features = {}
    for state in _FEATURES:
        all_features[state] = [_NEGATIVE_FLAG] * len(sequence)

    features = [_NEGATIVE_FLAG] * len(sequence)

    for i, seq in _bin(sequence, size=7):
        """
        Classify sequence in one of the disorder states

        See Fig 7 from Das & Pappu (2013)
        https://www.pnas.org/doi/full/10.1073/pnas.1304749110
        """
        f_plus = f_minus = 0

        for aa in seq:
            if aa in _POSITIVE_RESIDUES:
                f_plus += 1
            elif aa in _NEGATIVE_RESIDUES:
                f_minus += 1

        f_plus /= len(seq)
        f_minus /= len(seq)
        # fraction of charged residues
        fcr = f_plus + f_minus
        # net charge per residue
        ncpr = abs(f_plus - f_minus)

        state = None
        if fcr > 0.35:
            if ncpr <= 0.35:
                # Strong polyampholyte
                state = "Polyampholyte"
            elif f_plus > 0.35:
                # Strong positive polyelectrolyte
                state = "Positive Polyelectrolyte"
            elif f_minus > 0.35:
                # Strong negative polyelectrolyte
                state = "Negative Polyelectrolyte"
            else:
                raise ValueError(f"{seq}: {f_plus}, {f_minus}")
        else:
            # Weak polyampholyte/polyelectrolyte
            cysteine = proline = glycine = polar = 0
            for aa in seq:
                if aa == "C":
                    cysteine += 1
                elif aa == "P":
                    proline += 1
                elif aa == "G":
                    glycine += 1
                elif aa in {"N", "Q", "S", "T"}:
                    polar += 1

            if cysteine / len(seq) >= 0.32:
                state = "Cysteine-rich"
            elif proline / len(seq) >= 0.32:
                state = "Proline-rich"
            elif glycine / len(seq) >= 0.32:
                state = "Glycine-rich"
            elif seg_states[i]:
                state = "Low complexity"
            elif polar / len(seq) >= 0.32:
                state = "Polar"

        if state:
            all_features[state][i] = _POSITIVE_FLAG

    for state in reversed(_FEATURES):
        states = "".join(all_features[state])
        states = dilate(states, max_length=5)
        states = erode(states, max_length=5)

        # index = str(_FEATURES.index(state) + 1)
        for i, flag in enumerate(states):
            if flag == _POSITIVE_FLAG:
                features[i] = state

    return features


def _bin(sequence: str, size: int):
    """
    Split sequence in overlapping bins of size `size`
    :param sequence:
    :param size:
    :return:
    """
    n = math.floor((size - 1) / 2)
    seq_length = len(sequence)
    for i in range(seq_length):
        if i < n:
            yield i, sequence[n::-1] + sequence[1:n + 1]
        elif i + n < seq_length:
            yield i, sequence[i - n:i + n + 1]
        else:
            x = seq_length - (i + n) + 1
            yield i, sequence[i - n:] + sequence[::-1][1:x + 1]


def is_enriched(sequence: str, residues: set, threshold: float = 0.32):
    return len([aa for aa in sequence if aa in residues]) / len(sequence) >= threshold


def content_count(regions):
    return reduce(lambda x, t: x + (t[1] - t[0] + 1), regions, 0)
