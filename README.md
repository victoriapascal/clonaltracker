## Introduction

ClonalTracker is a tool to pairwise compare two vancomycin-resistant **Enterococcus faecium** (VRE) genomes 
that have been sequenced by Illumina-paired end sequencing. 
It compares the **van** gene type, the transposon type and the whole genome and makes a 
prediction whether the genomes are clonal or contain an identical transposon type or are unrelated.

## Installation

- Download ClonalTracker from gitlab (by using git clone or downloading a zip file): https://gitlab.com/victoriapascal/clonaltracker

```
git clone https://gitlab.com/victoriapascal/clonaltracker.git
```

- Install ClonalTracker dependencies by conda
```
conda env create -f clonaltracker/ct_env.yml ct_env
```

- Install all the dependencies one by one. First create a new conda environment:

```
conda create -n ct_env
```

```
conda activate ct_env
```

Then, install all the dependencies:

```
conda install poppunk==2.4.0
```

```
conda install pp-sketchlib==1.7.4
```

```
conda install -c bioconda isescan
```
Tested with ISEScan version v2.20

```
conda install -c bioconda ragtag
```
Tested with RagTag version v2.1.0

```
conda install -c bioconda clinker-py
```
Tested with clinker version v0.0.24

```
conda install -c bioconda mash
```
Tested with Mash version 2.3

```
conda install -c conda-forge -c bioconda bakta
```
For a full bakta installation, you also need to download its database by typing the following command:

```
bakta_db download --output <output-path> --type [light|full]

```
Tested with Bakta version 1.5.1 and full database. Please, check their repository for more information: https://github.com/oschwengers/bakta


## How to run ClonaTracker

- Activate your environment
```
conda activate ct_env
```

- Run ClonalTracker
	- Fasta files should be in the same folder
	- Path where to find both fasta files
```
python clonaltracker/clonaltracker.py /hpc/dla_mm/vpascalandreu/data/vanB_fastas/E7314.fasta /hpc/dla_mm/vpascalandreu/data/vanB_fastas/E7317.fasta  /path/to/bakta_db/
```

## Under-the-hood

### Input

ClonalTracker expects two assembled genomes as input (generated with e.g. SPAdes)

### Workflow

ClonalTracker runs the following tools:
- Blastn to detect the van type
- PopPUNK and MASH for whole genome comparison and contextualization in a larger dataset.
- Blastn, RagTag, ISEScan and Clinker to do transposon typing

Depending on the computational resources, ClonalTracker will take around 10 minutes to complete.

### Output

ClonalTracker outputs all the results and intermidiate files the different tools produce, which include interactive clinker HTML output, summary HTML file, log file, FASTA sequences...

## Test data

..

## Citation

