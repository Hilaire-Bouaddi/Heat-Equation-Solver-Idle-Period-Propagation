import matplotlib.pyplot as plt 
import numpy as np
import sys

with open(sys.argv[1]) as f:
    contents = f.read()
    NL = float(contents.splitlines()[0].split()[0])
    NH = float(contents.splitlines()[0].split()[1])
    N_ITER = float(contents.splitlines()[0].split()[2])

    T = []
    
    for n in [1]:
        print(contents.splitlines()[NH*n:NH*(n+1)])
        
