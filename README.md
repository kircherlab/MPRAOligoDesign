# Snakemake workflow: MPRAOligoDesign

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.20.0-brightgreen.svg)](https://snakemake.bitbucket.io)
![Tests](https://github.com/kircherlab/MPRAOligoDesign/workflows/Tests/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/mpraoligodesign/badge/?version=latest)](https://mpraoligodesign.readthedocs.io/en/latest/?badge=latest)

Workflow to generate an oligo design for MPRA experiments. It has multiple ways to generate sequences:

1. Variants (vcf) + regions (bed)
2. Variants (vcf)
3. Regions (bed)
4. Sequnces (fasta)

All 4 strategies can be combined later into one final design files. This allows you to put controls as well as your design in one workflow.

When sequences are generated it filters them with multiple check:
1. Homopolymer size (for all)
2. EcoRI or SbfI restriction site  (for all)
3. Simple repeats (not pure sequences)
4. TSS site overlap (not pure sequences)
5. CTCF overlap (not pure sequences)

## Authors

* Max Schubach (@visze), Berlin Institute of Health (BIH), [Computational Genome Biology](https://kircherlab.bihealth.org)

## Overview

![MPRA Oligo Design overview](docs/MPRAOligoDesign.png)

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above).

### Step 1: Obtain a copy of this workflow

[Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

### Step 2: Configure workflow

Configure the workflow according to your needs is via a config file. An example conmfig can be found in the `config/` folder. Adjust `config_test_samples.yml` to configure the workflow execution. For each input strategy you use (see above) You have to specify tab separated files with the (unique) sample id and where the input files are located in your file system. You can use relative or full paths (full paths recommended). relative means relative to your execution folder.

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

You can run the workflow in any directory. You just need read access to the cloned repository. This makes it easier to use the same workflow on multiple projects.

Test your configuration by performing a dry-run via

    snakemake --use-conda -n --configfile config.yaml --snakefile <path/tpo/cloned/workflow>/workflow/Snakefile

Execute the workflow locally via

    snakemake --use-conda --cores $N -configfile config.yaml --snakefile <path/tpo/cloned/workflow>/workflow/Snakefile

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100 -configfile config.yaml --snakefile <path/tpo/cloned/workflow>/workflow/Snakefile


If you not only want to fix the software stack but also the underlying OS, use

    snakemake --use-conda --use-singularity -configfile config.yaml --snakefile <path/tpo/cloned/workflow>/workflow/Snakefile

in combination with any of the modes above.
See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

### Step 5: Investigate results

After successful execution there should be a `results` folder. Under `results/final_design/` all sample ids are listed with the corresponding design file `results/final_design/<sample_id>/design.fa.gz`. A combined design file can be found here: `results/final_design/design.fa.gz`

If you want to get an overview why several sequences/regions are removed you can check the file `results/oligo_design/<sample_id>/design_<variants/sequences/regions>_filter.log`. `variants/sequences/regions` depends on the input you used.

### Step 6: Commit changes

Whenever you change something, don't forget to commit the changes back to your github copy of the repository:

    git commit -a
    git push

### Step 7: Obtain updates from upstream

Whenever you want to synchronize your workflow copy with new developments from upstream, do the following.

1. Once, register the upstream repository in your local copy: `git remote add -f upstream git@github.com:snakemake-workflows/MPRAOligoDesign.git` or `git remote add -f upstream https://github.com/snakemake-workflows/MPRAOligoDesign.git` if you do not have setup ssh keys.
2. Update the upstream version: `git fetch upstream`.
3. Create a diff with the current version: `git diff HEAD upstream/master workflow > upstream-changes.diff`.
4. Investigate the changes: `vim upstream-changes.diff`.
5. Apply the modified diff via: `git apply upstream-changes.diff`.
6. Carefully check whether you need to update the config files: `git diff HEAD upstream/master config`. If so, do it manually, and only where necessary, since you would otherwise likely overwrite your settings and samples.


### Step 8: Contribute back

In case you have also changed or added steps, please consider contributing them back to the original repository:

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the original repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to your local system, to a different place than where you ran your analysis.
3. Copy the modified files from your analysis to the clone of your fork, e.g., `cp -r workflow path/to/fork`. Make sure to **not** accidentally copy config file contents or sample sheets. Instead, manually update the example config files if necessary.
4. Commit and push your changes to your fork.
5. Create a [pull request](https://help.github.com/en/articles/creating-a-pull-request) against the original repository.

## Testing

Test cases are in the subfolder `.test`. They are automatically executed via continuous integration with [Github Actions](https://github.com/features/actions).



# Run oligo check scripts


## Create reference folder

mkdir reference

### Reference genome

```bash
cd reference
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz; bgzip hg38.fa
samtools faidx hg38.fa.gz
cd ..
```

### Simple repeat annotation

```bash
cd reference
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz
zcat simpleRepeat.txt.gz | cut -f 2- | sort -k1,1 -k2,2n -k3,3n | bgzip -c > simpleRepeat.bed.gz
tabix -p bed simpleRepeat.bed.gz
cd ..
```

### TSS positions

```bash
cd reference
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz
zcat gencode.v42.annotation.gtf.gz | awk 'BEGIN{ FS="\t"; OFS="\t" }{ if ($3 == "transcript") { if ($7 == "+") { print $1,$4-1,$4,".",0,"+" } else { print $1,$5-1,$5,".",0,"-" } } }' | sort -k1,1 -k2,2n -u | bgzip -c > TSS_pos.bed.gz
tabix -p bed TSS_pos.bed.gz
cd ..
```

## CTCF sites (ChIP peak + motif)

- Previously identified CTCF sites are in the data folder. Other definitions should be considered.

```bash
data/CTCF-MA0139-1_intCTCF_fp25.hg38.bed.gz
tabix -p bed data/CTCF-MA0139-1_intCTCF_fp25.hg38.bed.gz
```

## Applying available scripts

### Checking VCF files for minimal variant representation and correct reference allele

```bash
for i in *.vcf.gz; do zcat $i | python VCFrefFix.py --refCheck | bgzip -c > ${i/.vcf.gz/.fix.vcf.gz}; done
```

### Generating centered Ref/Alt sequences from VCF file

```bash
for i in *.vcf.fix.gz; do zcat $i | python varsToFrag.py | bgzip -c > ${i/.vcf.fix.gz/snvs.fa.gz}; done
```

### Design from pre-defined sequences with including barcodes in synthesis
- Note there is some deactivated tiling code based on BED files in this source code file
- Misses filter for TSS or CTCF, but probably easily added using similar code as for simple repeats
- Recently, we would no longer design with barcodes, but add barcodes through amplification; here only the outer sequences are added to the design 

```bash
./design.py -s input1.snvs.fa.gz,input2.snvs.fa.gz -r control_regions.fa.gz -o new_design.tsv 
```

