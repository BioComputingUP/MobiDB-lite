#!/usr/bin/env python3
"""
MobiDB-lite, 3.8.1, Jul 2020

By Marco Necci, Damiano Piovesan & Silvio C.E. Tosatto
BiocomputingUP lab, Padua, Italy

MobiDB-lite executes 8 different disorder predictors, collects the outputs and
calculates a consensus. The consensus is generated by measuring predictors
agreement. 5 out of 8 predictors must agree to assign disorder state to a
residue. Then a mathematical morphology (MM) dilation/erosion processing is
applied. Finally short regions are filtered out.

For further details on how to call mobidb_lite.py, call::

    mobidb_lite.py --help

For further details on requirements and troubleshooting, **see readme.md**
"""

import os
import json
import logging
import warnings
import configparser
from itertools import groupby

# relative imports
import mdblib.cli as cli
# import mdblib.plot as plot
import mdblib.logger as logger
from mdblib.protein import Protein
from mdblib.setdirs import set_pred_dir
from mdblib.streams import OutStream, InStream
from mdblib.consensus import MobidbLiteConsensus, SimpleConsensus
from mdblib.outformats import InterProFormat, ExtendedFormat, Mobidb3Format, CaidFormat

# Suppress warnings
warnings.filterwarnings('ignore')


class MobidbLite(object):
    """
    MobiDB-Lite application.

    MobiDB-Lite application launcher. Manages predictors execution,
    consensus computation and output formatting based on parameters.
    """
    outgroups = {'0': "main", '1': "main", '2': "mobidb3", '3': "caid"}

    def __init__(self, fasta, launchdir=None, conf=None, architecture='64', threads=0, outfile=None,
                 outfmt=0, skip_features=False, outmult_by='acc=', outmultsep=',',
                 parse_acc=False, force_consensus=False):

        cd = launchdir if launchdir else os.path.dirname(os.path.realpath(__file__))
        cd = os.path.abspath(cd)
        conf = self.read_config(conf if conf else os.path.join(cd, 'config.ini'))
        self.infname = fasta
        self.instream = None
        self.outstream = None

        # Set BINX directories
        self.bin_dirs = set_pred_dir(cd, conf.items('bin_directories'))
        self.active_preds = set(self.bin_dirs.keys())
        self.thresholds = {p: float(t) for p, t in dict(conf.items('thresholds')).items()}
        self.outgroup = self.outgroups[str(outfmt)]

        self.parse_acc = parse_acc
        self.skip_features = skip_features
        self.force_consensus = force_consensus
        self.outfmt = outfmt
        self.outmult_by = outmult_by
        self.outmultsep = outmultsep
        self.additional_data = None

        self.preds = self.run(fasta, architecture=architecture, threads=threads, outfile=outfile)

    def stream(self):
        """
        Stream MobidbLite output to outstream

        Consume generator produced by MobidbLite.run() and stream output to selected output-stream.
        Output-stream is selected by user at MobidbLite instantiation time (stdout by default).
        """
        for prot, s_cons, r_cons, m_cons in self.preds:
            outobj = self.fmt_output(prot.acc, prot.uniprot_acc, prot.seq, prot.preds,
                                     s_cons=s_cons, r_cons=r_cons, m_cons=m_cons)

            if outobj.isnone is False:
                self.outstream.write('{}\n'.format(outobj))

    @staticmethod
    def read_config(config):
        # Parse config file
        config_parser = configparser.ConfigParser()
        config_parser.optionxform = str
        config_parser.read(config)

        return config_parser

    def parse_json(self):
        for line in self.instream:
            doc = json.loads(line)
            acc, seq = doc["accession"], doc["sequence"]
            # remove field accession to prevent overwriting of existing accession field when injecting in output
            del doc["accession"]
            # del doc["sequence"]
            self.additional_data = doc
            yield acc, seq

    def parse_fasta(self):
        """
        Given a fasta file. yield tuples of header, sequence
        """
        faiter = (x[1] for x in groupby(self.instream, lambda line: line[0] == ">"))
        for header in faiter:
            # drop the ">"
            header = next(header)[1:].strip()
            # join all sequence lines to one.
            seq = "".join(s.strip() for s in next(faiter))
            logging.debug(header)
            yield header, seq

    def parse_input(self):
        """
        Sort out the correct parsing behaviour for input type
        :return: call to the parsing function
        :rtype: iterator
        """
        filename_elements = set(self.infname.split("."))
        fasta_extensions = {"fasta", "fa", "mfa", "mfasta", "multifasta", "multifa"}
        json_extensions = {"json", "mjson", "multijson", "multi-json"}
        if filename_elements & fasta_extensions:
            return self.parse_fasta()
        elif filename_elements & json_extensions:
            return self.parse_json()
        else:
            raise ValueError(
                "Unrecognized input-filename extension: '{}'. Expected one of {}".format(self.infname.split(".")[-1],
                                                                                         fasta_extensions | json_extensions))

    def fmt_output(self, acc, uacc, seq, preds, s_cons, r_cons, m_cons):
        output = None
        multi_acc = None

        # split fasta header (stored in acc) for multiple accessions
        if self.outmult_by is not None and self.outmult_by in acc:
            multi_acc = acc.split(self.outmult_by)[-1].split(self.outmultsep)

        # overwrite acc with a parsed accession if asked to
        if self.parse_acc:
            acc = uacc if uacc else acc

        # generate output based on selected output format
        if self.outfmt == 0:
            output = InterProFormat(acc, m_cons, _features=self.skip_features)

        elif self.outfmt == 1:
            output = ExtendedFormat(acc, m_cons, r_cons, _multi_accs=multi_acc)

        elif self.outfmt == 2:
            output = Mobidb3Format(acc, seq, m_cons, s_cons, preds, _multi_accs=multi_acc,
                                   injection=self.additional_data)

        elif self.outfmt == 3:
            output = CaidFormat(acc, seq, m_cons, preds, _multi_accs=multi_acc)

        return output

    def calc_consensus(self, predictions, sequence):
        simple_c = None
        relaxed_c = None
        mobidblite_c = MobidbLiteConsensus(predictions, sequence,
                                           pappu=True if self.outfmt == 2 else False,
                                           force=self.force_consensus)

        if self.outfmt == 1:
            relaxed_c = SimpleConsensus(predictions, sequence, force=self.force_consensus, threshold=.375)


        if self.outfmt == 2:
            simple_c = SimpleConsensus(predictions, sequence, force=self.force_consensus)



        return simple_c, relaxed_c, mobidblite_c

    def run(self, fasta, architecture, threads, outfile):

        logging.debug('outfmt: %i outgroup: %s', self.outfmt, self.outgroup)

        with InStream(fasta) as self.instream, OutStream(outfile) as self.outstream:
            # Parse input Fasta
            # for acc, sequence in self.parse_fasta():
            for acc, sequence in self.parse_input():
                # run predictors
                logging.debug('Current input %s', acc)

                with Protein(acc, sequence) as protein:
                    predictions = protein.run_predictors(
                        outgroup=self.outgroup,
                        active_preds=self.active_preds,
                        bin_dirs=self.bin_dirs,
                        thresholds=self.thresholds,
                        architecture=str(architecture),
                        processes=threads)

                    if predictions:
                        smp_consensus, rlx_consensus, mdbl_consensus = self.calc_consensus(predictions, sequence)
                        yield protein, smp_consensus, rlx_consensus, mdbl_consensus


if __name__ == "__main__":
    # Get dir where this piece of code is
    scriptdir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    # Parse command line arguments
    cli_args = cli.arg_parser(scriptdir)
    # Set logger
    logger.set_logger(cli_args.log, cli_args.logLevel)
    # Instantiate and run MobiDB-Lite application
    MobidbLite(cli_args.fastaFile,
               launchdir=scriptdir,
               conf=cli_args.conf,
               outfile=cli_args.outFile,
               outfmt=cli_args.outputFormat,
               architecture=cli_args.architecture,
               threads=cli_args.threads,
               skip_features=cli_args.skipFeatures,
               outmultsep=cli_args.multiplySeparator,
               outmult_by=cli_args.multiplyOutputBy,
               parse_acc=cli_args.parseAccession,
               force_consensus=cli_args.forceConsensus).stream()
