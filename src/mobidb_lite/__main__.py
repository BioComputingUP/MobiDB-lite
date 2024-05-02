import os
import sys
from argparse import ArgumentParser, FileType
from tempfile import gettempdir

from mobidb_lite.consensus import run


def main():
    script = os.path.relpath(__file__)

    description = "A consensus-based predictor of intrinsically " \
                  "disordered regions in proteins."
    parser = ArgumentParser(prog=f"python {os.path.basename(script)}",
                            description=description)
    parser.add_argument("infile", nargs="?", default="-",
                        type=FileType("rt", encoding="UTF-8"),
                        help="A file of sequences in FASTA format.")
    parser.add_argument("outfile", nargs="?", default="-",
                        type=FileType("wt", encoding="UTF-8"),
                        help="Write the output of infile to outfile.")
    parser.add_argument("--force", action="store_true", default=False,
                        help="Generate consensus as long as at least "
                             "one predictor did not fail.")
    parser.add_argument("--no-seg", dest="run_seg",
                        action="store_false", default=True,
                        help="Do not indentify domains "
                             "of low complexity with SEG.")
    parser.add_argument("--round", action="store_true", default=False,
                        help="Round scores before threshold checks, "
                             "like MobiDB-lite.")
    parser.add_argument("--tempdir", metavar="DIRECTORY", default=gettempdir(),
                        help=(f"Directory to use for temporary files, "
                              f"default: {gettempdir()}."))
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of parallel threads, default: 1.")
    args = parser.parse_args()

    root = os.path.abspath(os.path.dirname(script))
    bindir = os.path.join(root, "bin")

    with args.outfile as outfile:
        for seq_id, regions in run(args.infile, bindir, args.threads,
                                   force=args.force,
                                   round=args.round,
                                   seg=args.run_seg,
                                   tempdir=args.tempdir):
            if regions is None:
                sys.stderr.write(f"error in {seq_id}\n")
                continue

            for start, end, feature in regions:
                outfile.write(f"{seq_id}\t{start}\t{end}\t{feature}\n")


if __name__ == "__main__":
    main()
