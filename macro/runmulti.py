#!/usr/bin/env python3
# Run a program many times in parallel (small job pool).
#
# Each job is   target [args...] <job index>      (job index = 0 .. njob-1)
#
# Usage
#   python3 runmulti.py target title njob [ncore]
#   python3 runmulti.py -v target title njob -j 4 -- extra args
#
#   target : executable to run (e.g. ./test_rutherford_inspect)
#   title  : first argument given to the target (output file title)
#   njob   : number of jobs to run
#   ncore  : maximum number of jobs running at the same time (default 4)
#
# The script returns when every job is finished. Exit code is 0 if all the
# jobs succeeded and 1 if any of them failed (or 130 when interrupted with
# Ctrl-C, in which case all the running jobs are killed).

import argparse
import os
import signal
import subprocess
import sys
import time


def run_pool(cmd_base, njob, ncore, verbose=False, poll=0.2):
    """Run cmd_base + [str(i)] for i in range(njob), at most ncore at once.
    Returns the list of return codes indexed by job number."""
    codes = [None] * njob
    running = {}          # job index -> Popen
    next_job = 0
    try:
        while next_job < njob or running:
            while next_job < njob and len(running) < ncore:
                if verbose:
                    print("RUN  job %d (%d running)" % (next_job, len(running) + 1), flush=True)
                try:
                    running[next_job] = subprocess.Popen(cmd_base + [str(next_job)], start_new_session=True)
                except OSError as e:
                    print("FAILED to start job %d: %s" % (next_job, e), file=sys.stderr, flush=True)
                    codes[next_job] = 127
                next_job += 1
            for i in list(running):
                rc = running[i].poll()
                if rc is not None:
                    codes[i] = rc
                    del running[i]
                    if verbose or rc != 0:
                        print("DONE job %d (exit %d, %d running)" % (i, rc, len(running)), flush=True)
            if running:
                time.sleep(poll)
    except KeyboardInterrupt:
        print("KILL ALL PROCESS", file=sys.stderr, flush=True)
        for i, p in running.items():
            try:
                os.killpg(p.pid, signal.SIGKILL)   # the job's whole process group
            except OSError:
                p.kill()
        for i, p in running.items():
            p.wait()
            if verbose:
                print("JOB TERMINATED : %d" % i, flush=True)
        raise
    return codes


def main(argv=None):
    ap = argparse.ArgumentParser(description="Run 'target title <job index>' for njob jobs, ncore at a time.")
    ap.add_argument("target", help="executable to run")
    ap.add_argument("title", help="first argument of the target")
    ap.add_argument("njob", type=int, help="number of jobs")
    ap.add_argument("ncore", type=int, nargs="?", default=4, help="max simultaneous jobs (default 4)")
    ap.add_argument("-j", "--cores", type=int, dest="cores", help="same as ncore (overrides the positional one)")
    ap.add_argument("-q", "--quiet", action="store_true", help="do not print job start/finish messages")
    ap.add_argument("extra", nargs="*", help="extra arguments put between title and job index")
    a = ap.parse_args(argv)

    ncore = a.cores if a.cores is not None else a.ncore
    if a.njob < 0 or ncore < 1:
        ap.error("njob must be >= 0 and ncore >= 1")

    try:
        codes = run_pool([a.target, a.title] + a.extra, a.njob, ncore, not a.quiet)
    except KeyboardInterrupt:
        return 130
    failed = [i for i, c in enumerate(codes) if c != 0]
    if failed:
        print("%d/%d jobs FAILED: %s" % (len(failed), a.njob, failed), file=sys.stderr)
        return 1
    if not a.quiet:
        print("all %d jobs finished" % a.njob)
    return 0


if __name__ == "__main__":
    sys.exit(main())
