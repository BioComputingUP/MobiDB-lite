import logging

# relative imports
from mdblib.states import States


class Prediction(States):
    def __init__(self, method, scores, threshold, types=None,
                 include_in_mobidblite=False, invert=False):

        self.method = method
        self.scores = [round(x, 3) for x in (1-s if invert is True else s for s in scores)]
        self.threshold = threshold
        self.types = types
        self.include_in_mobidblite = include_in_mobidblite

        super(Prediction, self).__init__(self.scores_to_states())

    def __repr__(self):
        return "Prediction(method='{}', thr={}, scores=[{}..])".format(
            self.method, self.threshold, ' '.join(map(str, self.scores[:5])))

    def has_correct_length(self, seq, force=False):
        if (force is False and len(self.scores) == len(seq)) or force is True:
            return True
        else:
            logging.debug(
                'length difference | %s | len: %i pred: %s | %s',
                self.method, len(seq), len(self.scores), seq)
            return False

    def scores_to_states(self, tags=(1, 0)):
        """
        Threshold a continuous score into binary states

        :param tags: positive, negative labels for states
        :return: binary states
        """
        return [tags[0] if score >= self.threshold else tags[1] for score in self.scores]

    # def regions_to_set(self):
    #     """Get list of ID amino-acids positions from region list
    #
    #     :return: Unique amino-acid index of ID regions
    #     :rtype: set
    #     """
    #     positions = set()
    #     for start, end, _ in self.to_regions():
    #         positions.update(range(start, end + 1))
    #     return positions
    #
    # def regions_to_states(self, length, tags=('D', 'S'), reg_startindex=1):
    #     """Represent ID states as a string from a list ID regions
    #
    #     :param length: Length of the protein sequence (num of amino-acids)
    #     :type length: int
    #     :param tags: couple of values: positive match tag, negative match tag. Order counts
    #     :type tags: tuple
    #     :param reg_startindex: start index of the input regions (default: 1)
    #     :type reg_startindex: int
    #     :return: Order/Disorder states of the amino-acids of a protein
    #     :rtype: str
    #     """
    #     states = ""
    #     positions = [pos - reg_startindex for pos in self.regions_to_set()]
    #
    #     for i in range(0, length):
    #         if i in positions:
    #             states += tags[0]
    #         else:
    #             states += tags[1]
    #
    #     return states
