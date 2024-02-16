# CapSeq SNV calling pipeline 

Uses a snakemake pipeline to go from bams to mafs. Uses sage for SNV
calling. Includes some filtering steps as well. 

More details should be added here

Need to have snakemake env setup
```
mamba create -n snakemake snakemake
```


you need to copy over the ensemble reources into the cloned repo resources and unzip it
you can add it to gitignore, as the files are too big to push to github

```
Add 
cp /gsc/www/bcgsc.ca/downloads/morinlab/hmftools-references/ensembl_data_cache/38.zip cappseq_snv_pipeline/resources/ensemble_data_cache/
cd cappseq_snv_pipeline/resources/ensemble_data_cache/
unzip 38.zip

```