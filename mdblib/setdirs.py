import os
import logging


def set_pred_dir(cd, bindirs):
    src_dirs = dict()

    binxdir = os.path.join(cd, 'binx')
    packages_in_binx = os.listdir(binxdir)

    bindirs = dict(bindirs)

    missing = set(packages_in_binx) - set(bindirs.keys())
    if missing:
        logging.warning("%s %s been disabled", missing, 'has' if len(missing) == 1 else 'have')

    for pred, directory in bindirs.items():

        if not directory:
            default_bin = os.path.join(binxdir, pred)

            logging.debug(
                '{}: bin directory not set in config.ini; using default path: {}'.format(
                    pred.upper(), default_bin))

            src_dirs[pred] = default_bin
        else:
            src_dirs[pred] = os.path.abspath(directory)

    return src_dirs
