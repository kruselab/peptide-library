# peptide-library
This repository contains data and code described in Hurley *et al.* 2025 (currently unpublished).

Raw Illumina sequencing data, as well as processed files containing peptide sequences and their corresponding read counts, can be found in `data/`.

Code used to analyze NGS data can be found in `code/NGS.R`. A merged paired-end NGS dataset can be analyzed by running the function `run_NGS_pipeline("example_dataset.fastq.gz")`.

Code used to perform rarefaction-extrapolation analysis can be found in `code/rarefaction_analysis.R`. Multiple output CSV files from the NGS analysis pipeline above should be placed in the same directory, then this analysis can be performed using the function `run_rarefaction_analysis("example_directory/")`.

Code used to perform peptide library diversity and redundancy simulations can be found in `code/diversity_analysis.R`. As this process has fairly high computational and storage requirements for large libraries, clones are generated and counted in batches using a modified merge sort. Note also that this process was performed individually for each peptide length, after which uniques were summed to generate a final diversity estimate.
