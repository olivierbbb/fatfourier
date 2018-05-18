#!/usr/bin/env python3

import subprocess
import os
import sys

MIN_LENGTH = 20
MAX_LENGTH = 28

def iter_on_cmd(cmd):
    print('Running ' + cmd + '...')
    for i in range(MIN_LENGTH, MAX_LENGTH):
        full_cmd = cmd + " -l " + str(i)
        print(full_cmd)
        subprocess.run([full_cmd], shell=True)

def main():
    for file in os.scandir('bin/'):
        if not file.name.startswith('benchmark'):
            continue

        if file.name == 'benchmark_ff':
            for threads_count in [1, 2, 4, 8]:
                iter_on_cmd('bin/benchmark_ff -t {}'.format(threads_count))
        else:
            iter_on_cmd('bin/' + file.name)

main()
