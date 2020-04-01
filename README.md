# koji_mei_pipeline
https://github.com/shohei-kojima/koji_mei_pipeline


# prerequisites
- BLAST  2.8.1+ or later
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
-b /path/to/your_sample.bam \
-fa /path/to/hg38.fa \
-fai /path/to/hg38.fa.fai \
-fadb /path/to/hg38 \
-rep /path/to/humrepsub.fa \
-repout /path/to/hg38.fa.out \
-cov 35 \
-p 4


