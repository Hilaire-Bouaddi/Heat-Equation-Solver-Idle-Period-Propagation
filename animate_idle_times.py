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

            if i == 0 and j == 0:
                idle_times = np.zeros((N_ITER - 1, k, k))  

            values = list(map(int, contents.splitlines()[1].split()))
            
            for n in range(N_ITER - 1):
                idle_times[n][i][j] = values[n];

maximum_latency = np.max(np.max(np.max(idle_times, axis=2), axis=1), axis=0)
minimum_latency = np.min(np.min(np.min(idle_times, axis=2), axis=1), axis=0)

maximum_latency = 1e6

for n in range(N_ITER - 1): 
    plt.imshow(idle_times[n], cmap=plt.get_cmap('coolwarm'))#, vmin=minimum_latency, vmax=maximum_latency)
    plt.savefig("frames_idle_times/frame{}".format(n))

#fig, ax = plt.subplots()


"""
frames_names = []

for n in range(N_ITER):
    ax.cla()
    im = ax.imshow(T[n], cmap=plt.get_cmap('coolwarm'), vmin=0,vmax=100)
    ax.set_title("frame {}".format(n))
    # Note that using time.sleep does *not* work here!
    frame_name = "frames/img{}".format(n)
    plt.savefig(frame_name)
    frames_names.append(frame_name+".png")
    #plt.close()


fps = 1/dt
with imageio.get_writer('animation_parallel.gif', mode='I', fps=fps) as writer:
    for filename in frames_names:
        image = imageio.imread(filename)

"""
