import os
from subprocess import DEVNULL, PIPE, Popen


def run_anchor(bindir: str, file: str):
    cmd = [
        os.path.join(bindir, "anchor"),
        file
    ]
    proc = Popen(cmd, stdout=PIPE, stderr=DEVNULL, env={"ANCHOR_PATH": bindir},
                 encoding="utf-8")
    out, err = proc.communicate()

    if proc.returncode == 0:
        scores = []
        for line in out.splitlines():
            if line and line[0] != "#":
                pos, aa, score, state = line.split()
                scores.append(float(score))

        return scores

    return None


