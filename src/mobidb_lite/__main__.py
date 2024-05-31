import json
import os
import sys
from argparse import ArgumentParser, FileType
from tempfile import gettempdir

from mobidb_lite.consensus import run, content_count, _FEATURES, _MOBIDB_NAMES


def main():
    script = os.path.relpath(__file__)

    description = "A consensus-based predictor of intrinsically disordered regions in proteins."
    parser = ArgumentParser(prog=f"python {os.path.basename(script)}", description=description)
    parser.add_argument("infile", nargs="?", default="-", type=FileType("rt", encoding="UTF-8"),
                        help="A file of sequences in FASTA format.")
    parser.add_argument("outfile", nargs="?", default="-", type=FileType("wt", encoding="UTF-8"),
                        help="Write the output of infile to outfile.")
    parser.add_argument("--format", type=str, default='interpro',
                        help=f"Output format", choices=('interpro', 'mobidb')),
    parser.add_argument("--force", action="store_true", default=False,
                        help="Generate consensus as long as at least one predictor did not fail.")
    parser.add_argument("--no-seg", dest="run_seg", action="store_false", default=True,
                        help="Do not indentify domains of low complexity with SEG.")
    parser.add_argument("--nu", dest="calc_nu", action="store_true", default=False,
                        help="Calculate region compactness.")
    parser.add_argument("--round", action="store_true", default=False,
                        help="Round scores before threshold checks, like MobiDB-lite.")
    parser.add_argument("--tempdir", metavar="DIRECTORY", default=gettempdir(),
                        help=(f"Directory to use for temporary files, " f"default: {gettempdir()}."))
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of parallel threads, default: 1.")
    args = parser.parse_args()

    root = os.path.abspath(os.path.dirname(script))
    bindir = os.path.join(root, "bin")
    datadir = os.path.join(root, "data")

    with args.outfile as outfile:
        for seq_id, regions, scores in run(args.infile, bindir, datadir, args.threads,
                                   force=args.force,
                                   round=args.round,
                                   seg=args.run_seg,
                                   nu=args.calc_nu,
                                   tempdir=args.tempdir):
            if regions is None:
                sys.stderr.write(f"error in {seq_id}\n")
                continue

            seq_len = len(list(scores.values())[0])  # A random method

            if args.format == "interpro":
                for feature, region in regions.items():
                    if feature in _FEATURES:
                        for start, end in region:
                            outfile.write(f"{seq_id}\t{start}\t{end}\t{feature}\n")
                    elif feature == "mobidblite":
                        for start, end in region:
                            outfile.write(f"{seq_id}\t{start}\t{end}\t-\n")

            elif args.format == "mobidb":
                obj = {"acc": seq_id, "length": seq_len}
                # print(seq_id)
                for feature, region in regions.items():
                    cont_count = content_count(region)
                    cont_fraction = round(cont_count / seq_len, 3)
                    obj[_MOBIDB_NAMES[feature]] = {"regions": region, "content_count": cont_count, "content_fraction": cont_fraction}
                    if feature == "mobidblite":
                        obj[_MOBIDB_NAMES[feature]]["scores"] = scores[feature]
                outfile.write(json.dumps(obj) + "\n")


if __name__ == "__main__":
    main()
