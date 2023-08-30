# New RNA inference

This folder contains a computational pipeline that infers the fraction of new RNA per gene within a sample.

Describe how to design the yaml and run the pipeline.
Add a note on runtime and memory usage.


The inference script runs in sequential batches to keep memory usage low. The number of genes to be inferred per batch can be set in the Snakemake configuration file. Processing data using smaller batches requires less memory usage and does not significantly alter the runtime.

Average runtime will depend on the available hardware.

The pipeline has been benchmarked on a desktop computer with the following specs:
```
Hardware:
CPU: AMD Ryzen Threadripper 3960X, 24 core 48 threads
GPU: NVIDIA GeForce GTX 1080 Ti
Memory: 256 GB 

NVIDIA drivers:
Driver Version: 470.57.02
CUDA Version: 11.4 
```

Processing data from a large-scale experiment with 10 million data points, running on 30 parallel processes and performing 10.000 MCMC sinulations per batch finished in 96 hours and each batch required 130 GB RAM. 

## Usage
The pipeline consists of four separate steps:

### 1) Conversion counting from sequencing data
This step counts the number of observed mismatches to the reference genome observed in each sequencing reads and computes conversion and error probabilities. This is performed by the script `stranded_conversion_rates_and_probs.py`

#### Input
Configuration file in yaml-format.
Conversion counting starts from an aligned BAM-file. An example of how to generate one from raw sequencing data is found in `path_to_file`.

The supplied BAM file needs to contain read tags for genes and condition barcodes. By default, the pipeline looks for the gene tags `GE` and `GI` and barcode tag `BC`.

#### Output


#### 2) Combination and indexing of conversion count files
This step combines the output from the previous step

#### 3) Inference of $\pi_g$, the fraction of new RNA per gene in a treatment
#### 4) Collection of outdata and optional removal of intermediary folders

The resulting $\pi_g$ matrix contains the fraction of new RNA inferred per gene in each condition. Multiplying this 

Parameters for each step is entered into the `config.yaml` file. An example file can be found at `compound-screen/new_RNA_inference/config.yaml`. Once the yaml file is complete, place it in the desired run directory together with the Snakefile found at `compound-screen/new_RNA_inference/config.yaml`.

The full pipeline can then easily be started through: 

`snakemake --cores 1`

## Output

