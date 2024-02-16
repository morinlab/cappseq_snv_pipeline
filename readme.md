# CapSeq SNV calling pipeline 

Uses a snakemake pipeline to go from bams to mafs. Uses sage for SNV
calling. Includes some filtering steps as well. 

## Setup

Need to have snakemake env setup
```
mamba create -n snakemake snakemake
```


You need to copy over the ensemble resources into the cloned repo resources and unzip it
you can add it to gitignore, as the files are too big to push to github. Then make sure
the location is specified in the config, but it must be in your clone of this repo.

``` 
cp /gsc/www/bcgsc.ca/downloads/morinlab/hmftools-references/ensembl_data_cache/38.zip cappseq_snv_pipeline/resources/ensemble_data_cache/
cd cappseq_snv_pipeline/resources/ensemble_data_cache/
unzip 38.zip
```

## Troubleshooting missing samples

If you point the pipeline at a directory chockablock full of bams, and just give it a samplesheet and it cant find all the bams then some warnings will
be printed, and some files created to help you troubleshoot.


A way to weed out samples that the script cant find the bam files for you can import the SampleInfo object, and run it solo.
When it is run it will generate "missing_bams_log.txt" file that documents all the warnings about missing files. It also makes
two tsvs, one listing samples where bams could not be found, and another where normal bams could not be found. It assumes if there
is no UMI bam then there is no starting bam. 


```
import pandas as pd
import sys
sys.path.append('/home/kyakimovich/repos/cappseq_snv_pipeline') # update this path
import src.SampleInfo as SI
import datetime

TODAY = datetime.date.today().strftime("%Y-%m-%d")

# example paths to where bams are located
tumour_bam_dir = "/projects/rmorin/projects/DLBCL_montreal_thanos/PlasmaTwistCAPPSeq/BAMs/99-final/"
umi_bam_dir = "/projects/rmorin/projects/DLBCL_montreal_thanos/PlasmaTwistCAPPSeq/BAMs/04-umigrouped/"
normal_bam_dir = "/projects/rmorin/projects/DLBCL_montreal_thanos/NormalTargetted/targeted_seq_alignment_workflow/99-final/"

# spin up object, will print warnings
si_info = SI.SampleInfo(metadata_caps_nan,
              tumour_bam_dir,
              umi_bam_dir,
              normal_bam_dir)

# remove rows with missing samples 
miss_samples = pd.read_csv('missing_samples.tsv', sep='\t')
metadata_caps_nan_nomiss = metadata_caps_nan[~metadata_caps_nan['sample'].isin(miss_samples['missing_samples'])]

# write out new samplesheet without missing samples
samplesheet_nomissing.to_csv(f'samplesheet_nomissing_{TODAY}.tsv', sep='\t', index=False)

```
That or just figure out why your bam files are missing.