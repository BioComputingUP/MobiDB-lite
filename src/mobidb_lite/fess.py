import os
from subprocess import DEVNULL, PIPE, Popen


def parse(output):
    c = {"D", "O"}
    scores = []
    for line in output.splitlines():
        if line and line[0] in c:
            _, score = line.split("\t")
            scores.append(float(score))

    return scores


def run_fess(bindir: str, file: str):
    cmd = [
        os.path.join(bindir, "fess"),
        os.path.join(bindir, "model_fastss"),
        file,
        "/dev/null"
    ]
    proc = Popen(cmd, stdout=PIPE, stderr=DEVNULL, cwd=bindir, encoding="utf-8")
    out, err = proc.communicate()
    print(cmd)
    print(out)
    print(err)
    return parse(out) if proc.returncode == 0 else None

