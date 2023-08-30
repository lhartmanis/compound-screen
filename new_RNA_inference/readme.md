Describe how to design the yaml and run the pipeline.
Add a note on runtime and memory usage.


The inference script runs in sequantial batches to keep memory usage low. The number of genes to be inferred per batch can be set in the Snakemake configuration file. Processing data using smaller batches requires less memory usage and does not significantly alter the runtime.

Average runtime will depend on the available hardware.

The pipeline has been benchmarked on a desktop computer with the following specs:
```
CPU: AMD Ryzen Threadripper 3960X, 24 core 48 threads
GPU: NVIDIA GeForce GTX 1080 Ti
Memory: 256 GB 

NVIDIA drivers:
```

Processing data from a large-scale experiment with 10 million data points, running on 30 parallel processes and performing 10.000 MCMC sinulations per batch finished in 96 hours and each batch required 130 GB RAM. 
