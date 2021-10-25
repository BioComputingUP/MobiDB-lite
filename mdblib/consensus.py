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

    def calc_agreement(self, seq, threshold, ptype=None):
        """
        Compute agreement from a stack of predictions

        :param seq: amino acid sequence
        :type seq: str
        :param threshold: agreement threshold
        :type threshold: float
        :param ptype: prediction type to consider for agreement
        :type ptype: str
        """
        agreement = [0.0] * len(seq)
        included_predictors = 0

        for prediction in self.predictions_stack:

            if ptype is not None and prediction.has_correct_length(seq):
                if ptype in prediction.types:
                    logging.debug('%s | agreement: included', prediction.method)
                    logging.debug('%s | length: OK (%i)', prediction.method, len(seq))
                    included_predictors += 1
                    agreement = map(sum, zip(agreement, prediction.states))
            else:
                logging.warning('%s | agreement: excluded', prediction.method)

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
    def __init__(self, prediction_stack, seq, threshold=0.5):
        logging.debug('Generating Simple consensus')
        super(SimpleConsensus, self).__init__(prediction_stack)
        self.calc_agreement(seq, threshold, ptype='disorder')
        self.prediction.translate_states({1: 'D', 0: 'S'}, join_tr='')
        self.prediction.regions = self.prediction.to_regions(start_index=1, positivetag='D')


class MergeConsensus(Consensus):
    """
    Define a consensus merging all regions (e.g. for low complexity)
    """
    def __init__(self, prediction_stack, seq, threshold=0.1, ptype='disorder'):
        logging.debug('Generating Merge consensus')
        super(MergeConsensus, self).__init__(prediction_stack)
        self.calc_agreement(seq, threshold, ptype=ptype)
        self.prediction.regions = self.prediction.to_regions(start_index=1, positivetag=1)


class MobidbLiteConsensus(Consensus):
    """
    Define consensus featured by MobiDB-Lite as its prediction.
    """
    def __init__(self, prediction_stack, seq,
                 threshold=0.625, lencutoff=20, pappu=False, merge_features=True, keep_features_outside_idr=True):
        logging.debug('Generating MobiDB Lite consensus')
        super(MobidbLiteConsensus, self).__init__(prediction_stack)
        self.lencutoff = lencutoff
        self.seq = seq
        self.calc_agreement(seq, threshold, ptype='mobidblite')
        self.prediction.translate_states({1: 'D', 0: 'S'})
        self.prediction.math_morphology()
        self.prediction.merge_close_longidrs()
        self.prediction.apply_len_cutoff(self.lencutoff, tags=('D', 'S'), jn=True)
        self.prediction.regions = self.prediction.to_regions(start_index=1, positivetag='D',
                                                             len_thr=self.lencutoff)
        if self.prediction.regions:
            self.enriched_regions = self.get_region_features(merge=merge_features, only_in_idr=not keep_features_outside_idr)
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

    def get_region_features(self, window_size=9, feature_len_thr=10, merge=True, only_in_idr=True):
        """
        Look for sequence features within prediction.regions

        :return: prediction.regions extended with feature regions
        """

        def get_subregions_in_range(feature_iterable, start=0, end=None):
            subregions_in_range = []
            # transform features to regions in the disordered-region-range
            for feat_reg in States(feature_iterable[start: end]).to_regions(len_thr=feature_len_thr):
                if feat_reg[-1] != '0':
                    reg = [feat_reg[0] + start + 1, feat_reg[1] + start + 1, feature_desc[feat_reg[2]]]
                    subregions_in_range.append(reg)
            return subregions_in_range

        def append_subregions_to_enriched_regions(feature_iterable):
            # iterate over disordered regions
            for region in self.prediction.regions:
                start, end, _ = region
                # enriched_regions.append(region)
                if only_in_idr is True:
                    enriched_regions.extend(get_subregions_in_range(feature_iterable, start - 1, end))

            if only_in_idr is False:
                enriched_regions.extend(get_subregions_in_range(feature_iterable))

        if window_size % 2 == 0:
            logging.warning('Window size must be an odd integer instead it was passed: {}, using {} instead'.format(window_size, window_size+1))
            window_size += 1

        lc_regions = ''.join(map(str, next(p for p in self.predictions_stack if p.method == 'seg').states))

        features_raw = ['0'] * len(self.seq)
        features_merged = ['0'] * len(self.seq)
        enriched_regions = []

        seq = States(self.seq)
        for i, token in enumerate(seq.tokenize(n=window_size // 2 - 1)):
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

            elif len(lc_regions) == len(seq.states) and lc_regions[i] == '1':
                features_raw[i] = '7'

            elif token.is_enriched(['S', 'T', 'N', 'Q']) is True:
                features_raw[i] = '8'

        # slice forces copy, otherwise this assignment results in infinite loop
        enriched_regions = self.prediction.regions[:]

        # merge features hierarchically
        for feature_code in range(len(feature_desc), 0, -1):
            feature_code = str(feature_code)

            # apply math morph to single features
            f = States(features_raw)
            f.make_binary(active=feature_code)
            f.math_morphology(rmax=5, tags=(feature_code, '0'))

            if merge is True:
                # apply feature to seq positions
                for i, e in enumerate(f.states):
                    if e == feature_code:
                        features_merged[i] = feature_code
            else:
                # append features of each class without merging them
                append_subregions_to_enriched_regions(f.states)

        if merge is True:
            # append merged features
            append_subregions_to_enriched_regions(features_merged)

        return sorted(enriched_regions, key=lambda o: o[0])

