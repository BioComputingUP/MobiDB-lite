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
        self.skip_features = _features
        self.mdbl_consensus = _mdbl_consensus
        super(InterProFormat, self).__init__(_acc)

    def _get_output_obj(self):
        if self.mdbl_consensus.prediction.regions:
            if self.skip_features is False:
                out_obj = [[self.acc] + r for r in self.mdbl_consensus.enriched_regions]
            else:
                out_obj = [[self.acc] + r for r in self.mdbl_consensus.prediction.regions]

            self.isnone = False
            return [out_obj]

    def __repr__(self):
        if self.output:
            if self.skip_features is False:
                output = "\n".join(["{}\t{}\t{}{}".format(
                    ele[0], ele[1], ele[2], "\t{}".format(
                        ele[3]) if ele[3] != "D" else '') for ele in self.output[0]])
            else:
                output = "\n".join(["{}\t{}\t{}".format(
                        ele[0], ele[1], ele[2]) for ele in self.output[0]])
            return output
        else:
            return ""


class FastaFormat(Formatter):
    def __init__(self, _acc, _mdbl_consensus, _features=False, feature_desc=None):
        self.skip_features = _features
        self.mdbl_consensus = _mdbl_consensus
        self.feature_desc = {v: k for k, v in feature_desc.items()}
        super(FastaFormat, self).__init__(_acc)

    def _get_output_obj(self):
        if self.mdbl_consensus.prediction.regions:
            out_obj = {
                "accession": self.acc,
                "consensus": self.mdbl_consensus.prediction.states,
                "eregions": self.mdbl_consensus.enriched_regions if self.skip_features is False else None
            }
            self.isnone = False
            return [out_obj]

    def __repr__(self):
        if self.output:
            if self.skip_features is False:
                output = []
                for e in self.output:
                    c = list(e['consensus'])
                    for reg in e['eregions']:
                        if reg[2] != 'D':
                            c[reg[0]-1: reg[1]] = [self.feature_desc[reg[2]]] * (reg[1] - reg[0] + 1)
                    output.append(">{}\n{}".format(e['accession'], ''.join(c)))
                output = '\n'.join(output)
            else:
                output = '\n'.join([">{}\n{}".format(o['accession'], o['consensus']) for o in self.output])
            return output
        else:
            return ""


class Mobidb4Format(Formatter):

    def content_count(self, regions):
        return reduce(lambda x, t: x + (t[1] - t[0] + 1), regions, 0)

    def __init__(self, _acc, _seq, _mdbl_consensus,
                 _simple_consensus, _lowcomp_consensus, _single_predictions, **kwargs):
        self.seq = _seq
        self.seqlen = len(self.seq)
        self.mdbl_consensus = _mdbl_consensus
        self.simple_consensus = _simple_consensus
        self.single_predictions = _single_predictions
        self.lowcomplexity_consensus = _lowcomp_consensus
        self.injecting_data = kwargs.get("injection")
        super(Mobidb4Format, self).__init__(_acc, **kwargs)

        if self.multi_accessions:
            self.multiply_by_accession("accession")

    def _get_output_obj(self):
        out_obj = dict()

        if self.injecting_data is not None:
            out_obj.update(self.injecting_data)

        out_obj.setdefault("sequence", self.seq)

        # MobiDB-lite consensus and sub regions
        if self.mdbl_consensus.prediction.regions:
            regions = {}
            for r in self.mdbl_consensus.enriched_regions:
                r_type = "disorder" if r[2].startswith("D_") else r[2].replace("-", "_").replace(" ", "_").lower()
                regions.setdefault(r_type, []).append((r[0], r[1]))

            for r_type in regions:
                count = self.content_count(regions[r_type])
                obj_key = "prediction-{}-mobidb_lite{}".format(r_type, '_sub' if r_type != "disorder" else "")
                out_obj[obj_key] = {
                    'regions': regions[r_type],
                    'content_count': count,
                    'content_fraction': round(count / self.seqlen, 3)
                }
                if obj_key == "prediction-disorder-mobidb_lite":
                    out_obj[obj_key]["scores"] = self.mdbl_consensus.prediction.scores

        # Simple consensus
        if self.simple_consensus.prediction.regions:
            count = self.content_count(self.simple_consensus.prediction.regions)
            out_obj["prediction-disorder-th_50"] = {
                'regions': [(r[0], r[1]) for r in self.simple_consensus.prediction.regions],
                'content_count': count,
                'content_fraction': round(count / self.seqlen, 3)
            }

        if self.lowcomplexity_consensus.prediction.regions:
            count = self.content_count(self.lowcomplexity_consensus.prediction.regions)
            out_obj["prediction-low_complexity-merge"] = {
                'regions': [(r[0], r[1]) for r in self.lowcomplexity_consensus.prediction.regions],
                'content_count': count,
                'content_fraction': round(count / self.seqlen, 3)
            }

        # Single predictions
        for prediction in self.single_predictions:
            if prediction.valid is True:
                regions = [(r[0], r[1]) for r in prediction.to_regions(start_index=1, positivetag=1)]
                count = self.content_count(regions)

                if regions:
                    if 'disorder' in prediction.types:
                        out_obj["prediction-disorder-{}".format(prediction.method)] = {
                            'regions': regions,
                            'content_count': count,
                            'content_fraction': round(count / self.seqlen, 3)
                        }
                    elif 'lowcomp' in prediction.types:
                        out_obj["prediction-low_complexity-{}".format(prediction.method)] = {
                            'regions': regions,
                            'content_count': count,
                            'content_fraction': round(count / self.seqlen, 3)
                        }
                    elif 'bindsite' in prediction.types:
                        out_obj["prediction-lip-{}".format(prediction.method)] = {
                            'regions': regions,
                            'content_count': count,
                            'content_fraction': round(count / self.seqlen, 3)
                        }

                if 'rigidity' in prediction.types:
                    out_obj["prediction-rigidity-{}".format(prediction.method)] = {
                         'scores': prediction.scores
                    }

                if 'sspops' in prediction.types:
                    method, ptype = prediction.method.split('_')
                    out_obj["prediction-{}-{}".format(ptype, method)] = {
                        'scores': prediction.scores
                    }

        if out_obj:
            out_obj["length"] = self.seqlen

            if re.search("^UPI[A-F0-9]{10}$", self.acc):
                out_obj['uniparc'] = self.acc
            else:
                out_obj['acc'] = self.acc

            self.isnone = False

            return [out_obj]

    def __repr__(self):
        if self.output:
            return '\n'.join(json.dumps(oobj) for oobj in self.output)
        else:
            return ""


class CaidFormat(Formatter):
    def __init__(self, _acc, seq, _mdbl_consensus, _single_predictions, **kwargs):
        self.sequence = seq
        self.mdbl_consensus = _mdbl_consensus
        self.single_predictions = _single_predictions
        super(CaidFormat, self).__init__(_acc, **kwargs)

    def _get_output_obj(self):
        out_objs = list()
        if self.mdbl_consensus.prediction:
            self.mdbl_consensus.prediction.translate_states({'D': 1, 'S': 0})
            out_objs.append({
                "method": "mobidb_lite",
                "accession": self.acc,
                "sequence": self.sequence,
                "consensus": self.mdbl_consensus.prediction.states,
                "scores": self.mdbl_consensus.agreement
            })
            self.isnone = False
        for prediction in self.single_predictions:
            out_objs.append({
                "method": prediction.method,
                "accession": self.acc,
                "sequence": self.sequence,
                "consensus": prediction.states,
                "scores": prediction.scores
            })
            self.isnone = False
        return out_objs

    def __repr__(self):
        if self.output:
            return '\n'.join(['>{}\t{}\n{}'.format(
                o['accession'], o['method'],
                '\n'.join(['{}\t{}\t{}\t{}'.format(i + 1, *z) for i, z in enumerate(zip(o['sequence'], o['scores'], o['consensus']))])
            ) for o in self.output])

        else:
            return ""
