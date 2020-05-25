# koji_mei_pipeline
https://github.com/shohei-kojima/koji_mei_pipeline


# prerequisites
- BLAST  2.8.1+ or later (ncbi-blast-2.10.0+ cannot be used. Not compatible with biopython, 2020-05-25)
- bedtools  v2.26.0 or later
- Python  3.7 or later
- biopython  1.74 or later
- pysam  0.15.2 or later
- pybedtools  0.8.0 or later
- matplotlib  3.1.1 or later
- numpy  1.17.2 or later
- scipy  1.3.1 or later
- built-in python modules (os, sys, shutils, datetime, argsparse, glob, string, math, collections, itertools, multiprocessing, gzip, statistics)


# quick start with human samples
python main.py \
-overwrite -c ${cram_file} \
-fa GRCh38_full_analysis_set_plus_decoy_hla.fa \
-fadb GRCh38_full_analysis_set_plus_decoy_hla \
-mainchr ./lib/GRCh38DH_primary_plus_alt_ucsc_style.txt \
-rep ./lib/humrepsub.fa \
-repout hg38.fa.out \
-cov ${cov} -p 4


# requirements in a bam/cram file
SA:Z: tag
XA:Z: tag

# do not need in a bam/cram file (just a memo)
duplication tag
MC:Z: tag

