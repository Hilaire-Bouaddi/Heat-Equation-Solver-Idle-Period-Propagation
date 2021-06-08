# Heat Equation Solver & Idle Period Propagation

Final project for DD2356, Methods in High-Performance Computing

## The simulation 

### Building the code 

```cc -o parallel parallel.c```

### Running the code

```srun -n <n> ./parallel <L> <H> <h> <dt> <T_MAX> <k> [save_results] [save_idle_times]```

### Vizualization of the results as a GIF

After running with the command ```save_results```, run ```python animate.py results_parallel.py``` to generate a GIF.

## Monitoring of the Idle periods

After running with the command ```save_idle_times```, run ```python animate_idle_times.py <n>``` where n is the number of processes that generated the results to generate a GIF of heatmaps showing the latency at each node for every timestep.

