import os
import re
import base64
import codecs
import hashlib
import logging
import contextlib
from multiprocessing import Pool
from tempfile import NamedTemporaryFile

# relative imports
from mdblib import predictor


class Protein(object):
    """
    Defines a protein entity and manage operations relative to proteins.

    :param acc: identifier of the protein object
    :type acc: str
    :param seq: amino acid sequence of the protein object
    :type seq: str
    """
    uniprot_acc_pattern = re.compile(
        "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")

    def __init__(self, acc, seq):
        # passed attributes
        self.acc = acc
        self.seq = seq
        # computed attributes
        search = self.uniprot_acc_pattern.search(self.acc)
        self.uniprot_acc = search.group(0) if search else None
        self.secure_acc = self.acc.replace('|', '-').split()[0]
        # handle attributes
        self.seguid = None
        self.reprs = None
        self.preds = None

    def compute_seguid(self):
        """Generate protein sequence hash

        :return: hashed protein sequence
        :rtype: str
        """
        self.seq = codecs.latin_1_encode(self.seq.upper())[0]

        m = hashlib.sha1()
        m.update(self.seq)

        return base64.b64encode(m.digest()).rstrip("=")

    def __repr__(self):
        p = len(self.preds) if self.preds is not None else None
        return "Protein(acc='{}..', seq='{}..', predictions={})".format(self.acc[:10],
                                                                        self.seq[:10], p)

    def __enter__(self):
        """Generate temporary files representing a protein used as input of the predictors.

        Different predictors want different input formats. This function provides
        the 3 input format required from the protein sequence and accession. The 3
        formats required are:

        * disbin format::

            1
            sequence length
            sequence

        * flat format::

             sequence

        * fasta format::

            >accession
            sequence

        return: self
        """

        f_disbin = NamedTemporaryFile(
            delete=False, prefix="{}-disbin".format(self.secure_acc))
        f_flat = NamedTemporaryFile(
            delete=False, prefix="{}-flat".format(self.secure_acc))
        f_fasta = NamedTemporaryFile(
            delete=False, prefix="{}-fasta".format(self.secure_acc))

        f_disbin.write("1\n{}\n{}".format(len(self.seq), self.seq).encode('utf-8'))
        f_flat.write(self.seq.encode('utf-8'))
        f_fasta.write(">{}\n{}\n".format(self.acc, self.seq).encode('utf-8'))

        # set temporary files name
        logging.debug('Tempfiles generated')
        self.reprs = {
            'disbin': f_disbin.name,
            'flat': f_flat.name,
            'fasta': f_fasta.name
        }

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Delete temporary files.
        """
        with contextlib.suppress(FileNotFoundError):
            for fmt in self.reprs:
                os.remove(self.reprs[fmt])

    def run_predictors(self, outgroup, active_preds, bin_dirs, thresholds, architecture, processes):
        """Parallel call to predictors

        :param outgroup: tag indicating which predictors will be executed, can be configured in
            .ini file
        :type outgroup: str
        :param active_preds: predictors that have not been disabled by user
        :type active_preds: set
        :param bin_dirs: Directory of the predictor executables
        :type bin_dirs: dict
        :param thresholds: probability cutoff for discriminating ID from structure
            for each predictor
        :type thresholds: dict
        :param architecture: 32- or 64-bit OS architecture
        :type architecture: str
        :param processes: Number of worker processes of the process
            :py:class:`Pool` object
        :type processes: int
        """
        pool = Pool(processes) if processes > 0 else None

        preds = list()

        for subcl in predictor.Predictor.__subclasses__():
            if outgroup in subcl.groups and subcl.shared_name in active_preds:
                pred = subcl(self.reprs[subcl.intype], bin_dirs[subcl.shared_name],
                             architecture, thresholds)

                if pool is not None:
                    preds.append(pool.apply_async(pred.run))
                else:
                    logging.debug('Running predictor %s', pred.shared_name)
                    prediction = pred.run()
                    if prediction:
                        preds.extend(prediction)

        if pool is not None:
            pool.close()
            pool.join()

        if preds:
            if pool is not None:
                preds = self._unpack_pool_results(preds)

            self.preds = preds

        else:
            log_acc = self.uniprot_acc if self.uniprot_acc else self.secure_acc
            logging.error("%s | No predictors output", log_acc)

        return self.preds

    @staticmethod
    def _unpack_pool_results(pool_results):
        """Extract python data structures from pickled apply results

        :param pool_results: list of `multiprocessing.pool.ApplyResult`s
        :type pool_results: list

        :return: Unpacked apply results
        :rtype: dict
        """
        unpacked_results = list()
        for applyresult_obj in pool_results:

            try:
                if applyresult_obj.get():
                    for applyresult in applyresult_obj.get():
                        if applyresult:
                            unpacked_results.append(applyresult)
            except Exception as e:
                logging.warning(e)

        return unpacked_results
