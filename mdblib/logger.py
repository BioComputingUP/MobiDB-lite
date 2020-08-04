import logging
import sys


def set_logger(logfile, level):

    logging.basicConfig(level=logging.getLevelName(level),
                        format='%(asctime)s | %(module)-12s | %(levelname)-8s | %(message)s',
                        stream=open(logfile, "w") if logfile else sys.stderr)
