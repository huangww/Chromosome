"""
This is the script for chaning parameter and submiting job to the cluster
"""

import subprocess
import shlex
import fileinput
import os
import time
import numpy as np

def ready_to_submit():
    """ Inquire the queing system to determine the submit is ready or not
    :returns: Ture of False
    """
    inquire_process = subprocess.Popen(['qstat'], stdout=subprocess.PIPE)
    command = shlex.split("grep -w t")
    status = subprocess.Popen(command, stdin=inquire_process.stdout).communicate()
    return bool(status[0] is None)

def change_input(para, para_val, infile='input.in'):
    """ Change the input parameter

    :para: string, the parameter name to be changed
    :para_val: the parameter value to be set
    :infile: name of input file, default value is 'input.in'

    """
    for line in fileinput.input(infile, inplace=True):
        words = line.split()
        if len(words) == 3 and words[0] == para:
            line = line.replace(words[2], para_val)
        print line,
    return

def change_dt(dt_val, infile):
    """ Change the value of dt in the input file,
    also change other parameters correspondingly like tEnd

    :dt_val: the value to be set, string
    :infile: the name of input file, string

    """
    for line in fileinput.input(infile, inplace=True):
        words = line.split()
        if len(words) == 3 and words[0] == "dt":
            line = line.replace(words[2], str(dt_val))
        if len(words) == 3 and words[0] == "outputStep":
            line = line.replace(words[2], str(1e5))
        if len(words) == 3 and words[0] == "tEnd":
            line = line.replace(words[2], str(dt_val*1e4))
        print line,
    return

def find_dt(infile='input.in'):
    """ Find proper 'dt' to do the simulation

    :infile: input file name, string
    :returns: the value of dt, string

    """
    dt_upper = 1e-3
    dt_lower = 1e-6
    dt_val = dt_upper
    null_file = open(os.devnull, 'w')
    while dt_upper - dt_lower > 1e-6:
        change_dt(dt_val, infile)
        command = shlex.split("make run")
        try:
            exit_status = subprocess.call(command, stdout=null_file, stderr=subprocess.STDOUT)
            if exit_status == 0:
                dt_lower = dt_val
            else:
                dt_upper = dt_val
        except Exception as err:
            print err
            pass
        dt_val = (dt_upper + dt_lower) / 2.
    dt_string = ("%.1e" % dt_val).split('.')
    dt_val = dt_string[0] + dt_string[1][1:]
    return dt_val

def main():
    """ The main function
    :returns: TODO

    """
    start_time = time.time()
    para = "fExternal"
    para_range = np.linspace(0.1, 1, 0.1)
    for para_val in para_range:
        change_input(para, para_val)
        print para_val, find_dt('input.br.in')
    # print ready_to_submit()
    end_time = time.time()
    print "Total Running Time: ", end_time - start_time

if __name__ == "__main__":
    main()
