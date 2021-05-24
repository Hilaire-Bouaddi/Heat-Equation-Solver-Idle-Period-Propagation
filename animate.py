import matplotlib.pyplot as plt 
import numpy as np
import sys

with open(sys.argv[1]) as f:
    contents = f.read()
    NL = int(contents.splitlines()[0].split()[0])
    NH = int(contents.splitlines()[0].split()[1])
    N_ITER = int(contents.splitlines()[0].split()[2])

    T = []
    
    for n in range(N_ITER):
        T_n = []
        for line in contents.splitlines()[NH*n+1:NH*(n+1)+1]:
            T_n.append(list(map(float, line.split())))
        T.append(T_n)
    
T = np.array(T)

fig, ax = plt.subplots()

for n in range(N_ITER):
    ax.cla()
    im = ax.imshow(T[n], cmap=plt.get_cmap('coolwarm'), vmin=0,vmax=100)
    ax.set_title("frame {}".format(n))
    # Note that using time.sleep does *not* work here!
    plt.savefig("frames/frame {}".format(n))