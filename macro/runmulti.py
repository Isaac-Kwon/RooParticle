import sys, os, subprocess
from time import sleep

# Job Running Syntax
# target title(const char*) runnumber(Int)
#
# script running syntax
# python3 runmulti title njob ncore

class JobManager:
    def __init__(self, target, ncore=4, verbose=False):
        self.target  = target
        self.ncore   = ncore
        self.joblist = list()
        self.verbose = verbose
    def Run(self, args, njob=10):
        try:
            for jobN in range(njob):
                if len(self.joblist)<self.ncore:
                    if self.verbose:
                        print("NOW RUN : JOB NUMBER %d" %jobN)
                    self.joblist.append(subprocess.Popen([self.target]+args+[str(jobN)]))
                    sleep(1.0);
                if self.DeleteKilled():
                    pass
                else:
                    sleep(1.0);
        except KeyboardInterrupt:
            print("KILL ALL PROCESS")
            for i in range(len(self.joblist)):
                self.joblist[i].kill()
                self.joblist[i].wait()
                if self.verbose:
                    print("JOB TERMINATED : CURRENT JOB %d" %i)
            return
    def DeleteKilled(self):
        for i in range(len(self.joblist)):
            if self.joblist[i].poll() != None:
                self.joblist.pop(i)
                if self.verbose:
                    print("JOB KILLED : CURRENT JOB %d" %i)
                return True
        return False

if __name__ == "__main__":
    jm = JobManager(sys.argv[1], 10, True)
    jm.Run([sys.argv[2]], sys.argv[3])
