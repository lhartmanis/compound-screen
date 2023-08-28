# compound-screen
Contains code and analysis scripts for inferring fraction of new RNA from NASC-seq2 styled bulk experiments and analyzing gene expression patterns in new RNA. All expression data related to this project can be found at ArrayExpress with accession code **E-MTAB-13091**

### Dependencies
```
joblib (v1.0.1)
pysam (v0.17.0)
gtfparse(v1.2.1)
pandas(v1.1.4) 
numpy (v1.19.5) 
scipy (v1.6.3)
tensorflow (v2.6.0)
tensorflow_probability (v0.14.0)
```

### 1) Infer fraction new RNA
New RNA gets inferred per gene and condition using a Bayesian Markov chain Monte Carlo (MCMC) inference engine as described in *link to publication*.

This pipeline starts from a gene- and condition tagged BAM-file (as for example generated by [zUMIs](https://github.com/sdparekh/zUMIs)).
A Snakemake pipeline for inferring fraction of new RNA is provided under `new_RNA_inference`. The pipeline tries to speed up the computationally heavy MCMC inference by utilizing available GPUs. If no GPUs are found, all computations will run on the CPU. 

Code and examples of how to run the inference pipeline can be found in the `new_RNA_inference` folder

### 2) Analyze differentially expressed genes
