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

## How to run ClonaTracker

- Activate your environment
```
conda activate ct_env
```

- Run ClonalTracker
	- Fasta files should be in the same folder
	- Raw reads of both isolates should be in the same folder, and contain the *_R1.fastq.gz *_R2.fastq.gz suffix
```
python clonaltracker/clonaltracker.py /hpc/dla_mm/vpascalandreu/data/vanB_fastas/E7314.fasta /hpc/dla_mm/vpascalandreu/data/vanB_fastas/E7317.fasta /hpc/dla_mm/vpascalandreu/data/vanB_raw_reads_renamed2/
```

## Under-the-hood

### Input

ClonalTracker expects two paired-end fastq files and a corresponding assembly (generated with e.g. SPAdes)

### Workflow

ClonalTracker runs the following tools:
- Blastn to detect the van type
- TETyper, ISEScan, Clinker and BlastN for transposon typing
- PopPUNK and MASH for whole genome comparison and contextualization in a larger dataset.

Depending on the computational resources, ClonalTracker will take around 20 minutes to complete.

### Output

...

## Test data

..

## Command line options


Description of the command-line arguments 

## Citation

