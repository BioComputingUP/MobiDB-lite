import json
import gzip
import io
import logging
import os
import sys


class OutStream(object):
    def __init__(self, fn):
        self.filename = fn
        self.open_file = None

    def __enter__(self):
        # Set Output stream
        if self.filename:
            _, outfile_extension = os.path.splitext(self.filename)

            if outfile_extension == '.gz':
                self.open_file = gzip.open(self.filename, "wt")
            else:
                self.open_file = open(self.filename, "w")
        else:
            self.open_file = sys.stdout

        return self.open_file

    def __exit__(self, *args):
        if self.open_file != sys.stdout:
            self.open_file.close()

            if os.stat(self.filename).st_size == 0:
                os.remove(self.filename)


class InStream(object):
    def __init__(self, fn):
        self.filename = fn
        self.open_file = None

    def __enter__(self):
        if self.filename == '-':
            self.open_file = sys.stdin
        else:
            if os.path.isfile(self.filename):
                _, infile_extension = os.path.splitext(self.filename)
                if infile_extension == '.gz':
                    logging.debug('Opening Binary (gz \'ed) file, extension [%s]', infile_extension)
                    self.open_file = gzip.open(self.filename, 'rt')
                else:
                    logging.debug('Opening ASCII file, extension [%s]', infile_extension)
                    self.open_file = open(self.filename)
            else:
                if self.filename.lstrip()[0] == '>':
                    self.open_file = io.StringIO(self.filename.lstrip())
                else:
                    logging.critical("unrecognized input. Input is not a file NOR a fasta sequence")
                    raise ValueError("Unrecognized input")

        return self.open_file

    def __exit__(self, *args):
        if self.open_file != sys.stdin:
            self.open_file.close()
