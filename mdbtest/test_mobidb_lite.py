import json
import pytest
import numpy as np
# relative imports
from mobidb_lite import MobidbLite
from mdblib.prediction import Prediction, States
from mdblib.consensus import Consensus, SimpleConsensus, MobidbLiteConsensus
# data imports
from mdbtest.test_data import *


def test_prediction_hascorrectlength():
    scs = [0.1, 0.2, 0.5, 0.7, 0.9, 1]
    seq = 'a' * len(scs)
    assert Prediction('a', scs, 0.5).has_correct_length(seq)


@pytest.mark.parametrize("scores, thr, expected", [
    ([.1, .2, .5, .7, .9, 1], .5, [0, 0, 1, 1, 1, 1]),
    pytest.param([.1, .2, .4999, .7, .9, 1], .5, [0, 0, 0, 1, 1, 1], marks=[pytest.mark.xfail]),
    ([.1, .2, .499, .7, .9, 1], .5, [0, 0, 0, 1, 1, 1])
])
def test_prediction_states(scores, thr, expected):
    assert Prediction('a', scores, thr).states == expected


@pytest.mark.parametrize("states, trans, start, postag, lenthr, keepnone, expected", [
    ([0, 0, 0, 0, 1, 1, 1, 1], {}, 0, None, 1, False, [[0, 3, 0], [4, 7, 1]]),
    ([0, 0, 0, 0, 1, 1, 1, 1], {}, 1, None, 1, False, [[1, 4, 0], [5, 8, 1]]),
    ([0, 0, 0, 0, 1, 1, 1, 1], {0: 'S', 1: 'D'}, 1, None, 1, False, [[1, 4, 'S'], [5, 8, 'D']]),
    ([0, 0, 0, 0, 1, 1, 1, 1], {0: 'S', 1: 'D'}, 0, None, 1, False, [[0, 3, 'S'], [4, 7, 'D']]),
    ([0, 0, 0, 0, 1, 1, 1, 1], {}, 1, 1, 1, False, [[5, 8, 1]]),
    ([0, 0, 0, 0, 1, 1, 1, 1], {}, 1, 0, 1, False, [[1, 4, 0]]),
    ([0, 0, 0, 0, 1, 1, 1, 1, 1], {}, 1, None, 5, False, [[5, 9, 1]]),
    ([None, None, 0, 0, 1, 1, 1], {}, 1, None, 1, True, [[1, 2, None], [3, 4, 0], [5, 7, 1]]),
    ([None, None, 0, 0, 1, 1, 1], {}, 1, None, 1, False, [[3, 4, 0], [5, 7, 1]]),
])
def test_prediction_states_to_regions(states, trans, start, postag, lenthr, keepnone, expected):
    assert States(states).to_regions(keep_none=keepnone, start_index=start, translate_states=trans,
                                     positivetag=postag, len_thr=lenthr) == expected


@pytest.mark.parametrize("scores, dic, jn, expected", [
    ([0.1, 0.2, 0.5, 0.7, 0.9, 1], {0: 'S', 1: 'D'}, None, ['S', 'S', 'D', 'D', 'D', 'D']),
    ([0.1, 0.2, 0.5, 0.7, 0.9, 1], {0: 'S', 1: 'D'}, '', "SSDDDD"),
    ([0.1, 0.2, 0.5, 0.7, 0.9, 1], {0: 'S', 1: 'D'}, ',', "S,S,D,D,D,D"),
])
def test_translate_states(scores, dic, jn, expected):
    p = Prediction('a', scores, 0.5)
    p.translate_states(dic, join_tr=jn)

    assert p.states == expected


@pytest.mark.parametrize("scores, dic, jn, expected", [
    pytest.param([0.1, 0.2, 0.5, 0.7, 0.9, 1], {0: 'S', 1: 'D'}, 1, AttributeError),
    pytest.param([0.1, 0.2, 0.5, 0.7, 0.9, 1], {0: 'S'}, '', KeyError)
])
def test_fail_translate_states(scores, dic, jn, expected):
    with pytest.raises(expected):
        p = Prediction('a', scores, 0.5)
        p.translate_states(dic, join_tr=jn)


@pytest.mark.parametrize("tinput, expected", [
    ("SDSSSSDDSDDS", "SSSSSSDDDDDD"),
    ("DDDSSSDDDDDD", "DDDDDDDDDDDD"),
    ("SSSSSSDDSSSD", "SSSSSSSSSSSS")
])
def test_mathmorph(tinput, expected):
    s = States(tinput)
    s.math_morphology()

    assert s.states == expected


@pytest.mark.parametrize("tinput, rmax, expected", [
    ("SDSSSSDDSDDS", None, ValueError),
    ("SDSSSSDDSDDS", '3', ValueError),
])
def test_fail_mathmorph(tinput, rmax, expected):
    with pytest.raises(expected):
        s = States(tinput)
        s.math_morphology(rmax=rmax)


@pytest.mark.parametrize("tinput, expected", [
    ('D' * 21 + 'S' * 10 + 'D' * 21, 'D' * (21 + 10 + 21)),
    ('D' * 10 + 'S' * 10 + 'D' * 21, 'D' * 10 + 'S' * 10 + 'D' * 21),
    ('D' * 21 + 'S' * 11 + 'D' * 21, 'D' * 21 + 'S' * 11 + 'D' * 21),
])
def test_merge_close_longidrs(tinput, expected):
    s = States(tinput)
    s.merge_close_longidrs()

    assert s.states == expected


@pytest.mark.parametrize("seq, expected", [
    ('MAVMAPRTLLLLLSGALALTQTWAGSH', 'WC'),
    ('AAAAAAAAAAAAAAAAAAAAAAAAAAA', 'WC'),
    ('KKKKKKKKKKKKKKWIEQEGPEYWDQE', 'PA'),
    ('KKKKKKKKKKKKKKDDDDDDDDDDDDD', 'PA'),
    ('KKAAAAAAAAAAAAAAAAAAAAAADDD', 'WC'),
    ('DDDDDDDDDDDDAAAAAAAAAAAAAAA', 'NPE'),
    ('KKKKKKKKKKKKAAAAAAAAAAAAAAA', 'PPE'),
])
def test_get_disorder_class(seq, expected):
    assert States.get_disorder_class(seq) == expected


@pytest.mark.parametrize("x, y, pthr, cthr", [
    (8, 40, .5, .625),
    (80, 40, .5, .625),
    (800, 40, .5, .625),
    (800, 40, .2, .4),
    (800, 40, .9, .7),
    (800, 40, 1, 0),
    (800, 40, 0, 1),
    (800, 40, 0, 0),
    (800, 40, 1, 1),
])
def test_agreement(x, y, pthr, cthr):
    random_stack = np.random.rand(x, y).round(3)
    c = Consensus([Prediction('a', s, pthr, types=['a']) for s in random_stack.tolist()])
    c.calc_agreement('A'*y, cthr, ptype='a')
    assert c.prediction.states == ((random_stack >= pthr).mean(axis=0) >= cthr).astype(int).tolist()


@pytest.mark.parametrize("scores, cstates", [
    (escs, np.greater_equal(np.greater_equal(escs, .5).mean(axis=0), .5).astype(int))
])
def test_simple_consensus(scores, cstates):
    scores = np.array(scores)
    s = SimpleConsensus([Prediction('a', s, 0.5, types=['disorder']) for s in scores.tolist()],
                        'a' * scores.shape[1])
    s.prediction.translate_states({'D': 1, 'S': 0})
    assert s.prediction.states == cstates.tolist()


@pytest.mark.parametrize("scores, lcutoff, cstates", [
    (escs, 10, "DDDDDDDDDDSSSSDDDDDDDDDDDDDDSSSSSSSSSSSS"),
    (escs, 11, "SSSSSSSSSSSSSSDDDDDDDDDDDDDDSSSSSSSSSSSS"),
    (escs, 20, "SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS")
])
def test_mobidblite_consensus_states(scores, lcutoff, cstates):
    scores = np.array(scores)
    s = MobidbLiteConsensus(
        [Prediction('a', s, 0.5, types=['mobidblite']) for s in scores.tolist()],
        'a' * scores.shape[1], lencutoff=lcutoff)

    assert s.prediction.states == cstates


@pytest.mark.parametrize("scores, lcutoff, seq, cregs",  [
    (escs, 20, "MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPG", []),
    (escs, 11, "MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPG", [[15, 28, 'D_WC']]),
    (escs, 10, "MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPG", [[1, 10, 'D_WC'], [15, 28, 'D_WC']]),
    (escs, 10, "DDDDDDDDDDAAAAAKKKKKKKKKKKKKAAAAAAAAAAAA", [[1, 10, 'D_NPE'], [15, 28, 'D_PPE']]),
    (escs, 10, "DDDDAAAAAAAAAAAKKKKKKKDDDDDDAAAAAAAAAAAA", [[1, 10, 'D_NPE'], [15, 28, 'D_PA']]),
    (escs, 10, "DDDAAAAAAAAAAAAKKKKKKKDDDDDDAAAAAAAAAAAA", [[1, 10, 'D_WC'], [15, 28, 'D_PA']]),
])
def test_mobidblite_consensus_regions(scores, lcutoff, seq, cregs):
    sa = np.array(escs)
    ps = [Prediction('a', s, 0.5, types=['mobidblite']) for s in sa.tolist()]
    c = MobidbLiteConsensus(ps, seq, pappu=True, lencutoff=lcutoff)
    assert c.prediction.regions == cregs


@pytest.mark.parametrize("scores, lcutoff, seq, cregs",  [
    (ex1_scores, 20, ''.join(ex1_fasta.split('\n')[2:]), [[334, 365, 'D_WC'], [341, 359, 'Polar']]),
])
def test_region_features(scores, lcutoff, seq, cregs):
    ps = [Prediction('a', s, 0.2, types=['mobidblite']) for s in scores]
    c = MobidbLiteConsensus(ps, seq, pappu=True, lencutoff=lcutoff)

    assert c.enriched_regions == cregs

@pytest.mark.parametrize("states, n, expected",  [
    # small n
    ([1, 1, 0, 0, 1, 1, 1, 1], 1, [[1, 1, 1],
                                   [1, 1, 0],
                                   [1, 0, 0],
                                   [0, 0, 1],
                                   [0, 1, 1],
                                   [1, 1, 1],
                                   [1, 1, 1],
                                   [1, 1, 1]]),
    ([1, 1, 0, 0, 1, 1, 1, 1], 2, [[0, 1, 1, 1, 0],
                                   [1, 1, 1, 0, 0],
                                   [1, 1, 0, 0, 1],
                                   [1, 0, 0, 1, 1],
                                   [0, 0, 1, 1, 1],
                                   [0, 1, 1, 1, 1],
                                   [1, 1, 1, 1, 1],
                                   [1, 1, 1, 1, 1]]),
    ([1, 1, 0, 0, 1, 1, 1, 1], 3, [[0, 0, 1, 1, 1, 0, 0],
                                   [0, 1, 1, 1, 0, 0, 1],
                                   [1, 1, 1, 0, 0, 1, 1],
                                   [1, 1, 0, 0, 1, 1, 1],
                                   [1, 0, 0, 1, 1, 1, 1],
                                   [0, 0, 1, 1, 1, 1, 1],
                                   [0, 1, 1, 1, 1, 1, 1],
                                   [1, 1, 1, 1, 1, 1, 1]]),
    # n = len(states) - 1
    ([1, 1, 0, 0, 1, 1, 1, 1], 7, [[1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1],
                                   [1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1],
                                   [1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1],
                                   [1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1],
                                   [0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0],
                                   [0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0],
                                   [1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1],
                                   [1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1]]),
    # n > len(states) == len(states) - 1
    ([1, 1, 0, 0, 1, 1, 1, 1], 20, [[1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1],
                                    [1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1],
                                    [1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1],
                                    [1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1],
                                    [0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0],
                                    [0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0],
                                    [1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1],
                                    [1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1]]),
])
def test_tokenize(states, n, expected):
    s = States(states)
    assert list(s.tokenize(n=n)) == expected


@pytest.mark.parametrize("states, act, inct, expected",  [
    ('DDDDDSSDDDDCCC', 'D', 'S', 'DDDDDSSDDDDSSS'),
    ('DDDDDSSDDDDCCC', 'D', 'C', 'DDDDDCCDDDDCCC'),
    ('DDDDDSSDDDDCCC', 'A', 'S', 'SSSSSSSSSSSSSS'),
    (list('DDDDDSSDDDDCCC'), 'A', 'S', 'SSSSSSSSSSSSSS'),
])
def test_makebinary(states, act, inct, expected):
    s = States(states)
    s.make_binary(act, inct)
    assert s.states == expected


@pytest.mark.parametrize("fasta, outfmt, expected", [
    (ex1_fasta, 0, ex1_f0_stream),
    (ex1_fasta, 1, ex1_f1_stream),
    (ex1_fasta, 2, ex1_f2_stream),
    (ex1_fasta, 3, ex1_f3_stream),
])
def test_mobidblite_on_ex1(capsys, fasta, outfmt, expected):
    MobidbLite(fasta, outfmt=outfmt).stream()
    captured = capsys.readouterr()

    if outfmt == 0:
        assert captured.out == expected
    else:
        assert json.loads(captured.out) == json.loads(expected)
