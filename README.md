# MEGAnE
**MEGAnE**, **M**obile **E**lement search and **G**enotyping **An**alysis **E**nvironment, identifies and genotypes polymorphic mobile elements from short-read whole genome shotgun sequencing data (WGS). The current version does not support whole exome sequencing data nor is it tuned to detect somatic polymorphisms. The initial release of MEGAnE officially supports human and mouse datasets. However, we designed MEGAnE to allow analysis of other species, if the end user provides repeat library (e.g. consensus sequences from RepBase or Dfam).
  
MEGAnE (眼鏡 in Japanese) is pronunced like "mega" + "net." In Japanese, megane means a glass or glasses that fine-tunes our sight to see something or understand truth.  
  
**Currently MEGAnE is beta version. Please use this version at your own risk.**
  
# Citation
Mobile elements in human population-specific genome and phenotype divergence  
Shohei Kojima et al, bioRxiv 2022.03.25.485726; doi: https://doi.org/10.1101/2022.03.25.485726
  
# Installation
MEGAnE can be available as docker and Singularity containers from [dockerhub](https://hub.docker.com/r/shoheikojima/megane).  
We highly recommend to use such containers rather than preparing the required environment by yourself.  
  
```
# build for Singlarity
sudo singularity build MEGAnE.sif docker://shoheikojima/megane:v1.0.0.beta
# or 
singularity build --fakeroot MEGAnE.sif docker://shoheikojima/megane:v1.0.0.beta

# build for Docker
docker pull docker://shoheikojima/megane:v1.0.0.beta
```
  
# Input file
- MEGAnE can take a position-sorted BAM and CRAM file aligned by BWA-MEM and DRAGEN (in the case of DRAGEN, `-skip_unmapped` should be specified).  
- MEGAnE only supports paired-end WGS. Single-end WGS is not compatible. MEGAnE does not support WES.  
- We recommend to analyze WGS of 25x or higher depth, but it can also analyze 15x depth WGS by using the `-lowdep` option.  
- MEGAnE has a best performance with WGS of 150-bp or longer read length. We do not recommend to use WGS of less than 100-bp.  
- For more details, please see our [Wiki](https://github.com/shohei-kojima/MEGAnE/wiki) page.  
  
# In-depth usage
Please see our [Wiki](https://github.com/shohei-kojima/MEGAnE/wiki) page.  
  
# Quick usage for human WGS

### Step 0. Prepare MEGAnE k-mer file
- Before analyzing your BAM/CRAM files, you need to make MEGAnE k-mer files from your human reference genome (e.g. GRCh38DH, hs37d5, etc).  
- Usually, this takes ~10 min and requires ~50GB RAM.  
  
```
sif=/path/to/MEGAnE.sif

singularity exec ${sif} build_kmerset \
-fa /path/to/reference_human_genome.fa \
-prefix reference_human_genome \
-outdir megane_kmer_set
```
  
### Step 1. Call and genotype polymorphic MEs
- MEGAnE can analyze both BAM and CRAM format.  
- MEGAnE supports BAM/CRAM mapping to both GRCh37- and GRCh38-related genomes.  
- One 30x human WGS takes ~1 hour using 4 threads.  
  
```
sif=/path/to/MEGAnE.sif

# In the case of BAM file mapping to GRCh37-related genome (e.g. hs37d5, human_g1k_v37, hg19)
singularity exec ${sif} call_genotype_37 \
-i /path/to/input.bam \
-fa /path/to/reference_human_genome.fa \
-mk /path/to/megane_kmer_set/reference_human_genome.mk \
-outdir MEGAnE_result_test \
-sample_name test_sample \
-p 4

# In the case of CRAM file mapping to GRCh38-related genome (e.g. GRCh38DH, hg38)
singularity exec ${sif} call_genotype_38 \
-i /path/to/input.cram \
-fa /path/to/reference_human_genome.fa \
-mk /path/to/megane_kmer_set/reference_human_genome.mk \
-outdir MEGAnE_result_test \
-sample_name test_sample \
-p 4
```
  
### Step 2. Joint calling
- After the analysis of multiple BAM/CRAM files, you can make a joint call.  
- This will take several hours when merging 1000s samples.  
- MEGAnE supports joint calling from massive WGS (e.g. 10s of thousands). For more details, please see the [Wiki](https://github.com/shohei-kojima/MEGAnE/wiki) page.  
  
```
sif=/path/to/MEGAnE.sif

# first, list up samples (output directories from Step 1) you are going to merge
ls -d /path/to/[all_output_directories] > dirlist.txt

# merge non-reference ME insertions
singularity exec ${sif} joint_calling_hs \
-merge_mei \
-f dirlist.txt \
-fa /path/to/reference_human_genome.fa \
-cohort_name test

# merge reference ME polymorphisms
singularity exec ${sif} joint_calling_hs \
-merge_absent_me \
-f dirlist.txt \
-fa /path/to/reference_human_genome.fa \
-cohort_name test
```
  
### (Optional) Step 3. Make a joint call for haplotype phasing
- The step 2 above will generate two VCF files. This step merges the two files generated in the step 2.  
- When merging the two VCF files, MEGAnE removes multi-allelic variants.  
- This is particularly useful when you do haplotype phasing using MEGAnE's results.  
  
```
sif=/path/to/MEGAnE.sif

singularity exec ${sif} reshape_vcf \
-i /path/to/jointcall_out/[cohort_name]_MEI_jointcall.vcf \
-a /path/to/jointcall_out/[cohort_name]_MEA_jointcall.vcf \
-cohort_name test
```
  
### (Optional) Example haplotype phasing of MEGAnE's result  
- Here is one example of how to phase MEGAnE's result.  
- This step requires external softwares.  
  
```
threads=6
memory=30000

# first, generate a SNP VCF
snp_vcf=/path/to/SNP.vcf chr=1
exclude=/path/to/vcf_for_phasing/cohort_name_biallelic.bed.gz
chr=1

plink2 \
--threads ${threads} \
--memory ${memory} \
--vcf ${snp_vcf} \
--make-pgen \
--max-alleles 2 \
--mac 2 \
--hwe 1e-6 \
--chr ${chr} \
--indiv-sort natural \
--exclude bed0 ${exclude} \
--export vcf bgz \
--out SNP


# next, generate a ME VCF
me_vcf=/path/to/vcf_for_phasing/[cohort_name]_biallelic.vcf.gz

plink2 \
--threads ${threads} \
--memory ${memory} \
--vcf ${me_vcf} \
--make-pgen \
--max-alleles 2 \
--mac 2 \
--hwe 1e-6 \
--chr ${chr} \
--indiv-sort natural \
--vcf-half-call missing \
--export vcf bgz \
--out ME


# merge SNP and ME VCFs
cat SNP.vcf.gz > _SNP_ME.vcf.gz
zcat ME.vcf.gz | grep -v '#' | pigz -c >> _SNP_ME.vcf.gz
zcat _SNP_ME.vcf.gz | pigz -c > SNP_ME.vcf.gz
rm  _SNP_ME.vcf.gz

plink2 \
--threads ${threads} \
--memory ${memory} \
--vcf SNP_ME.vcf.gz \
--make-pgen \
--sort-vars \
--vcf-half-call missing \
--out SNP_ME

plink2 \
--threads ${threads} \
--memory ${memory} \
--pfile SNP_ME \
--export vcf bgz \
--out SNP_ME_sorted

tabix SNP_ME_sorted.vcf.gz


# haplotype-phasing
map=/path/to/genetic_map

shapeit4 \
--input SNP_ME_sorted.vcf.gz \
--map ${map} \
--region ${chr} \
--thread ${threads} \
--output SNP_ME.phased.vcf.gz
```

