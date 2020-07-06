import re
import json
import copy
from functools import reduce


class Formatter(object):
    def __init__(self, _acc, _multi_accs=None, **kwargs):
        self.acc = _acc.split()[0]
        self.multi_accessions = _multi_accs
        self.isnone = True
        self.output = self._get_output_obj()

    def multiply_by_accession(self, key):
        for i, acc in enumerate(self.multi_accessions):
            if i == 0:
                self.output[0][key] = acc
            else:
                copy_out = copy.copy(self.output[0])
                copy_out[key] = acc
                self.output.append(copy_out)

    def _get_output_obj(self):
        # handle
        pass


class InterProFormat(Formatter):
    def __init__(self, _acc, _mdbl_consensus, _features=False):
        self.features = _features
        self.mdbl_consensus = _mdbl_consensus
        super(InterProFormat, self).__init__(_acc)

    def _get_output_obj(self):
        if self.mdbl_consensus.prediction.regions:
            if self.features is False:
                out_obj = [[self.acc] + r for r in self.mdbl_consensus.enriched_regions]
            else:
                out_obj = [[self.acc] + r for r in self.mdbl_consensus.prediction.regions]

            self.isnone = False
            return [out_obj]

    def __repr__(self):
        if self.output:
            if self.features is False:
                output = "\n".join(["{}\t{}\t{}{}".format(
                    ele[0], ele[1], ele[2], "\t{}".format(
                        ele[3]) if ele[3] != "D" else '') for ele in self.output[0]])
            else:
                output = "\n".join(["{}\t{}\t{}".format(
                        ele[0], ele[1], ele[2]) for ele in self.output[0]])
            return output
        else:
            return ""


class ExtendedFormat(Formatter):
    def __init__(self, _acc, _mdbl_consensus, rlx_consensus, **kwargs):
        self.mdbl_consensus = _mdbl_consensus
        self.rlx_consensus = rlx_consensus
        super(ExtendedFormat, self).__init__(_acc, **kwargs)

        if self.multi_accessions:
            self.multiply_by_accession("accession")

    def _get_output_obj(self):
        if self.mdbl_consensus.prediction.regions:
            out_obj = {
                "accession": self.acc,
                "consensus": self.mdbl_consensus.prediction.states,
                "rlx_consensus": self.rlx_consensus.prediction.states,
                "regions": self.mdbl_consensus.prediction.regions
            }
            self.isnone = False
            return [out_obj]

    def __repr__(self):
        if self.output:
            return '\n'.join(json.dumps(oobj) for oobj in self.output)
        else:
            return ""


class Mobidb3Format(Formatter):
    def __init__(self, _acc, _seq, _mdbl_consensus,
                 _simple_consensus, _single_predictions, **kwargs):
        self.seq = _seq
        self.seqlen = len(self.seq)
        self.mdbl_consensus = _mdbl_consensus
        self.simple_consensus = _simple_consensus
        self.single_predictions = _single_predictions
        self.injecting_data = kwargs.get("injection")
        super(Mobidb3Format, self).__init__(_acc, **kwargs)

        if self.multi_accessions:
            self.multiply_by_accession("accession")

    def _get_output_obj(self):
        out_obj = dict()

        if self.injecting_data is not None:
            out_obj.update(self.injecting_data)

        out_obj.setdefault("sequence", self.seq)

        # MobiDB-lite consensus
        out_obj \
            .setdefault('mobidb_consensus', dict()) \
            .setdefault('disorder', dict()) \
            .setdefault('predictors', list()) \
            .append(
            {'method': 'mobidb_lite',
             'regions': self.mdbl_consensus.prediction.regions,
             'scores': self.mdbl_consensus.prediction.scores,
             'dc': reduce(
                 lambda x, t:
                 x + (t[1] - t[0] + 1),
                 self.mdbl_consensus.prediction.regions, 0.0) / self.seqlen
                 if self.mdbl_consensus.prediction.regions else 0.0})

        # MobiDB-lite consensus sub regions
        if self.mdbl_consensus.prediction.regions:
            out_obj \
                .setdefault('mobidb_consensus', dict()) \
                .setdefault('disorder', dict()) \
                .setdefault('predictors', list()) \
                .append(
                {'method': 'mobidb_lite_sub',
                 'regions': self.mdbl_consensus.enriched_regions_tags
                 })

        # Simple consensus
        out_obj.setdefault('mobidb_consensus', dict()) \
            .setdefault('disorder', dict()) \
            .setdefault('predictors', list()) \
            .append(
            {'method': 'simple',
             'regions': self.simple_consensus.prediction.regions,
             'dc': reduce(
                 lambda x, t:
                 x + (t[1] - t[0] + 1),
                 self.simple_consensus.prediction.regions, 0.0) / self.seqlen
                 if self.simple_consensus.prediction.regions else 0.0})

        # Single predictions
        for prediction in self.single_predictions:
            if any(t in prediction.types for t in ['disorder', 'lowcomp']):
                prediction.translate_states({1: 'D', 0: 'S'})
                out_obj \
                    .setdefault('mobidb_data', dict()) \
                    .setdefault('disorder', dict()) \
                    .setdefault('predictors', list()) \
                    .append(
                    {'method': prediction.method,
                     'regions': prediction.to_regions(start_index=1, positivetag='D')})

            if 'sspops' in prediction.types:
                method, ptype = prediction.method.split('_')

                out_obj \
                    .setdefault('mobidb_data', dict()) \
                    .setdefault('ss_populations', dict()) \
                    .setdefault('predictors', list()) \
                    .append(
                    {'method': method,
                     'type': ptype,
                     'scores': prediction.scores})

            if 'bindsite' in prediction.types:
                prediction.translate_states({1: 'D', 0: 'S'})

                out_obj \
                    .setdefault('mobidb_data', dict()) \
                    .setdefault('lips', dict()) \
                    .setdefault('predictors', list()) \
                    .append(
                    {'method': prediction.method,
                     'regions': prediction.to_regions(start_index=1, positivetag='D')})

        if out_obj:
            out_obj["length"] = self.seqlen

            if re.search("^UPI[A-F0-9]{10}$", self.acc):
                out_obj['uniparc'] = self.acc

            else:
                out_obj['accession'] = self.acc

            self.isnone = False

            return [out_obj]

    def __repr__(self):
        if self.output:
            return '\n'.join(json.dumps(oobj) for oobj in self.output)
        else:
            return ""


class CaidFormat(Formatter):
    def __init__(self, _acc, seq, _mdbl_consensus, _single_predictions, **kwargs):
        self.seq = seq
        self.mdbl_consensus = _mdbl_consensus
        self.single_predictions = _single_predictions
        super(CaidFormat, self).__init__(_acc, **kwargs)

        if self.multi_accessions:
            self.multiply_by_accession("accession")

    def _get_output_obj(self):
        out_obj = dict()

        out_obj \
            .setdefault('predictions', list()) \
            .append(
            {'method': 'mobidb_lite',
             'regions': self.mdbl_consensus.prediction.regions,
             'scores': self.mdbl_consensus.prediction.scores})

        for prediction in self.single_predictions:
            if 'disorder' in prediction.types:
                prediction.translate_states({1: 'D', 0: 'S'})

                out_obj \
                    .setdefault('predictions', list()) \
                    .append(
                    {'method': prediction.method,
                     'regions': prediction.to_regions(start_index=1, positivetag='D'),
                     'scores': prediction.scores})

            if 'sspops' in prediction.types:
                method, _ = prediction.method.split('_')

                out_obj \
                    .setdefault('predictions', list()) \
                    .append(
                    {'method': method,
                     'regions': prediction.to_regions(start_index=1,
                                                      positivetag=1,
                                                      translate_states={1: 'D', 0: 'S'}),
                     'scores': prediction.scores})

            if 'bindsite' in prediction.types:
                prediction.translate_states({1: 'D', 0: 'S'})

                out_obj \
                    .setdefault('predictions', list()) \
                    .append(
                    {'method': prediction.method,
                     'regions': prediction.to_regions(start_index=1, positivetag='D'),
                     'scores': prediction.scores})

        # if out_obj:
        out_obj['accession'] = self.acc
        out_obj['sequence'] = self.seq

        self.isnone = False
        return [out_obj]

    def __repr__(self):
        if self.output:
            return '\n'.join(json.dumps(oobj) for oobj in self.output)
        else:
            return ""
