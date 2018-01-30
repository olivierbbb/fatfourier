#!/usr/bin/env python3

import subprocess
import os
import sys

if len(sys.argv) != 3:
    print("usage: " + sys.argv[0] + " <bin/fft_to_test> <path/to/audio.wav>", file = sys.stderr)
    sys.exit()

fft_cmd = sys.argv[1]
path_to_audio = sys.argv[2]

#i = 20
i = 2
while i < 28:
    print("\n2^" + str(i))
    full_cmd = fft_cmd + " " + path_to_audio + " " +str(2**i)
    print(full_cmd)
    subprocess.run([full_cmd], shell=True)
    i += 1
