# About MEGAnE

# Installation
MEGAnE can be available as Docker and Singularity containers. We highly recommend to use such containers.  
  
```
# build for Singlarity 3
sudo singularity build MEGAnE_v1.0.0.sif docker://shoheikojima/megane:v1.0.0

# build for Docker
docker pull docker://shoheikojima/megane:v1.0.0
```
  
# In-depth usage
Please see our wiki page (under construction) or the instruction PDF file found at `docs/MEGAnE_v1.0_instruction.pdf`.  
  
# Quick usage for human WGS

### Step 0. Prepare MEGAnE k-mer file
- Before analyzing your BAM/CRAM files, you need to make MEGAnE k-mer files from your human reference genome (e.g. GRCh38DH, hs37d5, etc).  
- Usually, this takes ~10 min and requires ~50GB RAM.  
  
```
sif=/path/to/MEGAnE_[version].sif

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
sif=/path/to/MEGAnE_[version].sif

# In the case of BAM file mapping to GRCh37-related genome (e.g. hs37d5, human_g1k_v37, hg19)
singularity exec ${sif} call_genotype_37 \
-b /path/to/input.bam \
-fa /path/to/reference_human_genome.fa \
-mk /path/to/megane_kmer_set/reference_human_genome.mk \
-outdir test_run \
-sample_name test_sample \
-p 4

# In the case of CRAM file mapping to GRCh38-related genome (e.g. GRCh38DH, hg38)
singularity exec ${sif} call_genotype_38 \
-c /path/to/input.cram \
-fa /path/to/reference_human_genome.fa \
-mk /path/to/megane_kmer_set/reference_human_genome.mk \
-outdir MEGAnE_result_sample_1 \
-sample_name sample_1 \
-p 4
```
  
### Step 2. Joint calling
- After the analysis of multiple BAM/CRAM files, you can make a joint call.  
- This will take several hours when merging 1000s samples.  
  
```
sif=/path/to/MEGAnE_[version].sif

# first, list up samples you are going to merge
ls -d /path/to/MEGAnE_result_* > dirlist.txt

# merge non-reference ME insertions
singularity exec ${sif} build_kmerset \
-merge_mei \
-f dirlist.txt \
-fa /path/to/reference_human_genome.fa \
-cohort_name test

# merge reference ME polymorphisms
singularity exec ${sif} build_kmerset \
-merge_absent_me \
-f dirlist.txt \
-fa /path/to/reference_human_genome.fa \
-cohort_name test
```
  
### (Optional) Step 3. Make a joint call for imputation
- The step 2 above will generate two VCF files. This step merges the two files generated in the step 2.  
- When merging the two VCF files, MEGAnE removes multi-allelic variants.  
- This is particularly useful when you do haplotype-phasing using MEGAnE's results.  
  
```
sif=/path/to/MEGAnE_[version].sif

singularity exec ${sif} reshape_vcf \
-i /path/to/jointcall_out/[cohort_name]_MEI_jointcall.vcf \
-a /path/to/jointcall_out/[cohort_name]_MEA_jointcall.vcf \
-cohort_name test
```
  
### Example haplotype-phasing for MEGAnE's result  
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

