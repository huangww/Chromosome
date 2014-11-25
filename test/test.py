import os

files = os.listdir('.')

for fname in fles:
    os.rename(fname, fname + '.modified')


