"""
This is the script for chaning parameter and submiting job to the cluster
"""

import subprocess
import shlex
import fileinput
import time
import re
import numpy as np

def ready_to_submit():
    """ Inquire the queing system to determine the submit is ready or not
    :returns: Ture of False
    """
    inquire_process = subprocess.Popen(['qstat'], stdout=subprocess.PIPE)
    output = inquire_process.stdout.read()
    # status = ('qw' in output) or ('t' in output)
    status = bool(re.search(r'\bqw\b', output)) or bool(re.search(r'\bt\b', output))
    return not status


def job_submit(command):
    """ Sumit job to the queuing system

    :command: the command used to submit job, "qsub"

    """
    subprocess.call(shlex.split(command))
    return


def change_input(para, para_val, infile='input.in'):
    """ Change the input parameter

    :para: string, the parameter name to be changed
    :para_val: the parameter value to be set
    :infile: name of input file, default value is 'input.in'

    """
    for line in fileinput.input(infile, inplace=True):
        words = line.split()
        if len(words) == 3 and words[0] == para:
            line = line.replace(words[2], str(para_val))
        print line,
    return


def main():
    """ The main function

    """
    start_time = time.time()
    para = 'fExternal'
    para_range = np.linspace(0.1, 1.0, 10)
    for para_val in para_range:
        while not ready_to_submit():
            time.sleep(60)
        change_input(para, para_val)
        command = "qsub -t 1-100 -N stretch_coil_N100 job.sh"
        job_submit(command)
        print para_val
    end_time = time.time()
    print "Total Running Time: ", end_time - start_time

if __name__ == "__main__":
    main()
    # print ready_to_submit()
