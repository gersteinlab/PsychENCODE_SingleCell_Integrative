# PsychENCODE_snRNA_pipeline
The snRNA-seq processing pipeline was constructed based on published best practices and new benchmarking metrics for existing methods. The pipeline is mostly implemented in Python, except for the final cell-type annotation steps that involve the R-based program Azimuth. The overall workflow can be summarized in three main steps:
1. Count matrix generation, demultiplexing, and ambient RNA clean-up
2. Per-fastq set/sample processing using [Pegasus](https://pegasus.readthedocs.io/en/stable/)
3. Per-study aggregation of processed sample and cell-type annotation

Documentation is provided in the file ```Detailed_PsychENCODE_snRNAseq_pipeline.pdf```.

To clone the conda environments appropriate for running different parts of the pipeline, please use the ```Pegasus_spec-file.txt``` and ```R_spec-file.txt``` files to clone them (see instructions [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#building-identical-conda-environments)).
