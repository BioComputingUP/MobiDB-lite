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

