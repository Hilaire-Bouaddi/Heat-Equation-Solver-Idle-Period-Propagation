import matplotlib.pyplot as plt 
import numpy as np
import sys

N_THREADS = int(sys.argv[1])
k = np.sqrt(N_THREADS).astype(np.int64)

for i in range(k):
    for j in range(k):
        filename = "idle_times/coords_{},{}.txt".format(i, j)
        with open(filename) as f:
            contents = f.read();
            N_ITER = int(contents.splitlines()[0].split()[2])
            dt = float(contents.splitlines()[0].split()[3])

            if i == 0 and j == 0:
                idle_times = np.zeros((N_ITER - 1, k, k))  

            values = list(map(int, contents.splitlines()[1].split()))
            
            for n in range(N_ITER - 1):
                idle_times[n][i][j] = values[n];


for i in range(k): 
    for j in range(k):
        plt.plot(np.arange(0, dt * (N_ITER - 1), dt), idle_times[:, i, j]/np.max(idle_times[:, i, j]), linestyle="--")

latency_sum = np.sum(np.sum(idle_times, axis=1), axis=1)

plt.plot(np.arange(0, dt * (N_ITER - 1), dt), latency_sum[:]/np.max(latency_sum),linewidth=4, label="latencies summed")

plt.xlabel("time (s)")
plt.ylabel("latency")
plt.title("Normalized latency evolution per process")
plt.legend()


plt.savefig("plot_idle_times")

