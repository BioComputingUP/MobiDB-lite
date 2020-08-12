"""Contains classes used to build different type of consensus."""

import logging

# relative imports
from mdblib.prediction import Prediction
from mdblib.states import States


pappu_codes = {'PA': '1',
               'PPE': '2',
               'NPE': '3'}

feature_codes = {'C': '4',
                 'P': '5',
                 'G': '6',
                 'STNQ': '8'}

feature_desc = {'1': 'Polyampholyte',  # PA
                '2': 'Positive Polyelectrolyte',  # PPE
                '3': 'Negative Polyelectrolyte',  # NPE
                '4': 'Cystein-rich',  # CR
                '5': 'Proline-rich',  # PR
                '6': 'Glycine-rich',  # GR
                '7': 'Low complexity',  # LC
                '8': 'Polar'}  # PO

feature_tag = {'Polyampholyte': 'PA',  # PA
               'Positive Polyelectrolyte': 'PPE',  # PPE
               'Negative Polyelectrolyte': 'NPE',  # NPE
               'Cystein-rich': 'CR',  # CR
               'Proline-rich': 'PR',  # PR
               'Polar': 'PO',  # PO
               'Glycine-rich': 'GR',  # GR
               'Low complexity': 'LC'}  # LC


class Consensus(object):
    """
    General consensus class. Can calculate the raw agreement of predictors.
    """
    tag = 'mdb'

    def __init__(self, prediction_stack):
        self.predictions_stack = prediction_stack
        self.summed_states = None
        self.agreement = None
        self.prediction = None

    def __repr__(self):
        return "Consensus(methods='{}', {})".format(
            ' '.join(p.method for p in self.predictions_stack), self.prediction)

    def calc_agreement(self, seq, threshold, ptype=None, force_consensus=False):
        """
        Compute agreement from a stack of predictions

        :param seq: amino acid sequence
        :type seq: str
        :param threshold: agreement threshold
        :type threshold: float
        :param ptype: prediction type to consider for agreement
        :type ptype: str
        :param force_consensus: if True consensus computation is computed despite
            single predictors errors
        :type force_consensus: bool
        """
        agreement = [0.0] * len(seq)
        included_predictors = 0

        for prediction in self.predictions_stack:

            if ptype is not None and ptype in prediction.types:
                logging.debug('%s | agreement: included', prediction.method)

                if prediction.has_correct_length(seq, force_consensus):
                    logging.debug('%s | length: OK (%i)', prediction.method, len(seq))
                    included_predictors += 1
                    agreement = map(sum, zip(agreement, prediction.states))
            else:
                logging.debug('%s | agreement: excluded', prediction.method)


        if included_predictors != 0:
            self.summed_states = agreement
            agreement = [summed_states / included_predictors for summed_states in agreement]
            self.agreement = agreement
        else:
            logging.error('No predictions where included in agreement calculation')

        self.prediction = Prediction(self.tag, agreement, threshold)


class SimpleConsensus(Consensus):
    """
    Define a simple consensus based on an agreement threshold (default 0.5).
    """
    def __init__(self, prediction_stack, seq, threshold=0.5, force=False):
        logging.debug('Generating Simple consensus')
        super(SimpleConsensus, self).__init__(prediction_stack)
        self.calc_agreement(seq, threshold, ptype='disorder', force_consensus=force)
        self.prediction.translate_states({1: 'D', 0: 'S'}, join_tr='')
        self.prediction.regions = self.prediction.to_regions(start_index=1, positivetag='D')


class MergeConsensus(Consensus):
    """
    Define a consensus merging all regions (e.g. for low complexity)
    """
    def __init__(self, prediction_stack, seq, threshold=0.1, ptype='disorder', force=True):
        logging.debug('Generating Simple consensus')
        super(MergeConsensus, self).__init__(prediction_stack)
        self.calc_agreement(seq, threshold, ptype=ptype, force_consensus=force)
        self.prediction.regions = self.prediction.to_regions(start_index=1, positivetag=1)


class MobidbLiteConsensus(Consensus):
    """
    Define consensus featured by MobiDB-Lite as its prediction.
    """
    def __init__(self, prediction_stack, seq,
                 threshold=0.625, lencutoff=20, pappu=False, force=False):
        logging.debug('Generating MobiDB Lite consensus')
        super(MobidbLiteConsensus, self).__init__(prediction_stack)
        self.lencutoff = lencutoff
        self.seq = seq
        self.calc_agreement(seq, threshold, ptype='mobidblite', force_consensus=force)
        self.prediction.translate_states({1: 'D', 0: 'S'})
        self.prediction.math_morphology()
        self.prediction.merge_close_longidrs()
        self.prediction.apply_len_cutoff(self.lencutoff, tags=('D', 'S'), jn=True)
        self.prediction.regions = self.prediction.to_regions(start_index=1, positivetag='D',
                                                             len_thr=self.lencutoff)
        if self.prediction.regions:
            self.enriched_regions = self.get_region_features()
            self.enriched_regions_tags = [(r[0], r[1], 'D_' + feature_tag[r[2]]) for r in
                                          filter(lambda r: r[-1][:2] != 'D', self.enriched_regions)]
            if pappu is True:
                self.set_pappu_classes_per_region()

    def set_pappu_classes_per_region(self, reg_startindex=1):
        """
        Transform the status of regions appending the Pappu class. (setter)

        :param reg_startindex: start index of regions
        """
        for region in self.prediction.regions:
            start, end, status = region
            region_sequence = self.seq[start - reg_startindex: end - reg_startindex + 1]
            pappu_class = self.prediction.get_disorder_class(region_sequence)
            region[-1] = '{}_{}'.format(status, pappu_class)

    def get_region_features(self):
        """
        Look for sequence features within prediction.regions

        :return: prediction.regions extended with feature regions
        """

        lc_regions = ''.join(map(str, next(p for p in self.predictions_stack if p.method == 'seg').states))

        features_raw = ['0'] * len(self.seq)
        features_final = ['0'] * len(self.seq)
        enriched_regions = []

        seq = States(self.seq)
        for i, token in enumerate(seq.tokenize(n=7)):
            token = States(token)

            pappu_class = token.get_disorder_class(token.states)

            # assign features hierarchically
            if pappu_class != 'WC':
                features_raw[i] = pappu_codes[pappu_class]

            elif token.is_enriched(['C']) is True:
                features_raw[i] = '4'

            elif token.is_enriched(['P']) is True:
                features_raw[i] = '5'

            elif token.is_enriched(['G']) is True:
                features_raw[i] = '6'

            elif lc_regions[i] == '1':
                features_raw[i] = '7'

            elif token.is_enriched(['S', 'T', 'N', 'Q']) is True:
                features_raw[i] = '8'

            logging.debug(features_raw)

        # merge features hierarchically
        for feature_code in range(len(feature_desc), 0, -1):
            feature_code = str(feature_code)

            # apply math morph to single features
            f = States(features_raw)
            f.make_binary(active=feature_code)
            f.math_morphology(rmax=5, tags=(feature_code, '0'))

            # apply feature to seq positions
            for i, e in enumerate(f.states):
                if e == feature_code:
                    features_final[i] = feature_code

        for region in self.prediction.regions:
            start, end, _ = region
            enriched_regions.append(region)

            for feat_reg in States(features_final[start - 1: end]).to_regions(len_thr=15):
                if feat_reg[-1] != '0':
                    reg = [feat_reg[0] + start, feat_reg[1] + start, feature_desc[feat_reg[2]]]
                    enriched_regions.append(reg)

        return enriched_regions
