# Classes for processes, implementing some of the methods of
# multiprocessing.Process

from __future__ import with_statement
import os, time
import util

class Process(object):
    """A dummy class to be used when no parallel processes are
    possible/desired."""
    def __init__(self, target, args):
        self.target = target
        self.args = args
        self.exitcode = None

    def start(self):
        self.target(*self.args) 
        self.exitcode = 0

    def is_alive(self):
        """This process executes in serial, so if you get to call
        is_alive, the answer must be no."""
        return False

class SGEprocess(Process):
    """A class for shell commands executed using the Sun Grid Engine.
    Partly based on sge.py by Trevor Strohman
    http://ciir.cs.umass.edu/~strohman/code/sge.py
    """

    verbose = False

    def __init__(self, command, name, directory, priority=1, nslots=1):
        self.exitcode = None
        self.command = command
        self.name = name
        self.priority = priority
        self.nslots = nslots
        self.dir = directory

        self.outdir = os.path.join(self.dir, 'stdout')
        self.errdir = os.path.join(self.dir, 'stderr')
        self.scriptdir = os.path.join(self.dir, 'scripts')
        dirs = [self.dir, self.outdir, self.errdir, self.scriptdir]
        map(os.mkdir, filter(lambda(d): not os.path.exists(d), dirs))

        self.scriptfile = os.path.join(self.scriptdir, self.name)
        self.outfile = os.path.join(self.outdir, self.name)
        self.errfile = os.path.join(self.errdir, self.name)
        with open(self.scriptfile, 'w') as f:
            f.write(self.shell_script())

    def shell_script(self):
        """Construct the qsub submission shell script"""
        script = []
        script += ["#$ -N %s" % self.name]
        script += ["#$ -o %s" % self.outfile]
        script += ["#$ -e %s" % self.errfile]
        script += ["#$ -cwd"]
        script += ["#$ -V"]
        script += ["#$ -pe level%d.pe %d" % (self.priority, self.nslots)]
        script += [self.command]
        return '\n'.join(script) + '\n'
    
    def start(self):
        if self.verbose: print(self.command)
        os.system("qsub %s > /dev/null" % self.scriptfile)

    def is_alive(self):
        code = os.system( "qstat -j %s >& /dev/null" % (self.name) )
        if code != 0:
            self.set_exitcode()
            return False 
        return True

    def wait(self):
        interval = 5
        # I worry that there might be a delay between qsub and
        # appearance of job in qstat, so
        time.sleep(interval)
        while self.is_alive():
            time.sleep(interval)
            interval = min( 2 * interval, 60 )
            
    def set_exitcode(self):
        while not os.path.exists(self.errfile):
            if self.verbose: print('%s: waiting for error output' % self.name)
            time.sleep(10)
        self.exitcode = util.count_lines(self.errfile)

    def execute_in_serial(self):
        self.start()
        self.wait()
        if self.exitcode != 0:
            raise ProcessError('SGE job %s produced error output: %s' % (self.name, self.errfile))


class ProcessError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
