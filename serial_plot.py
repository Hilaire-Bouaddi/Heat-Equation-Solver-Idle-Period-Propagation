import numpy as np
import matplotlib.pyplot as plt

N_ITER_NH_NL = [
        600000,
        2400000,
        60000000,
        106134000,
        166333500,
        240000000
        ]

T = [
        [0.012808, 0.012205, 0.012056, 0.012001, 0.013174],
        [0.050498, 0.050691, 0.050074, 0.057313, 0.057728],
        [1.340528, 1.634543, 1.298753, 1.309342, 1.296952],
        [2.201734, 2.223768, 2.700633, 2.723344, 2.200620],
        [3.184282, 3.164972, 3.184333, 3.763918, 3.171684],
        [4.460560, 4.455190, 4.454526, 4.941094, 4.453768]
    ]

T = np.array(T) 

plt.errorbar(N_ITER_NH_NL, np.mean(T, axis=1), yerr=np.std(T, axis=1), marker='o')
plt.title("Execution time of the serial version of the code")
plt.xlabel("N_ITER * NH * NL")
plt.ylabel("time (s)")

plt.savefig("serial_plot")


