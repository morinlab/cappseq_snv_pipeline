#!/usr/bin/env snakemake

import os
import pandas as pd
import pyfaidx
import numpy
import snakemake

snakemake.utils.min_version("7")

import sys
sys.path.append(config['cappseq_snv_pipeline']['pipeline_dir'])
import src.SampleInfo as SI

# Load input files
samplesheet = pd.read_csv(config["cappseq_snv_pipeline"]["samplesheet"], sep="\t")

SampleInfo = SI.SampleInfo(samplesheet, 
                config["cappseq_snv_pipeline"]["tumour_bam_dir"],
                config["cappseq_snv_pipeline"]["umi_bam_dir"],
                config["cappseq_snv_pipeline"]["normal_bam_dir"],)
# can remove this by using object through out, do one day
samples = SampleInfo.samples
samplepaths = SampleInfo.samplepaths
normals = SampleInfo.normals
normal_ids = SampleInfo.normal_ids
sample_uncons = SampleInfo.sample_uncons


# Run variant calling
rule run_sage:
    input:
        tbam = lambda w: samplepaths[w.sample],
        nbam = lambda w: normals[w.sample],
    output:
        vcf = config["cappseq_snv_pipeline"]["base_dir"] + "/01-SAGE/{sample}/{sample}.sage.vcf"
    params:
        ref_genome = config["cappseq_snv_pipeline"]["ref_genome"],
        ref_genome_version = "38" if config["cappseq_snv_pipeline"]["ref_genome_ver"] == "GRCh38" else "37",  # Should make more robust
        # Panel regions
        hotspots_vcf = os.path.join(config["cappseq_snv_pipeline"]["pipeline_dir"], "resources/KnownHotspots.vcf.gz"),
        panel_regions = config["cappseq_snv_pipeline"]["capture_space"],
        ensembl = os.path.join(config["cappseq_snv_pipeline"]["pipeline_dir"], "resources/ensembl_cache/38/"),
        # Miscellaneous
        normal_name = lambda w: normal_ids[w.sample],
        max_depth = config["cappseq_snv_pipeline"]["max_depth"],
        min_map = config["cappseq_snv_pipeline"]["min_map_qual"],
        hard_vaf_cutoff = config["cappseq_snv_pipeline"]["tumor_min_vaf"],
        soft_vaf_cutoff = config["cappseq_snv_pipeline"]["tumor_soft_min_vaf"],
        min_norm_depth = 7
    threads: 4
    log:
        config["cappseq_snv_pipeline"]["base_dir"] + "/logs/{sample}.sage_run.log"
    conda:
        "envs/sage.yaml"
    shell:
        """SAGE -tumor {wildcards.sample} -tumor_bam {input.tbam} \
        -out {output.vcf} -ref_genome {params.ref_genome} \
        -ref_genome_version {params.ref_genome_version} \
        -hotspots {params.hotspots_vcf} -panel_bed {params.panel_regions} \
        -high_confidence_bed {params.panel_regions} \
        -ensembl_data_dir {params.ensembl} \
        -reference_bam {input.nbam} -reference {params.normal_name} \
        -max_read_depth {params.max_depth} -min_map_quality {params.min_map} \
        -hard_min_tumor_vaf {params.hard_vaf_cutoff} -hard_min_tumor_raw_alt_support 4 \
        -panel_min_tumor_vaf {params.soft_vaf_cutoff} -panel_min_germline_depth {params.min_norm_depth} \
        -hotspot_min_germline_depth {params.min_norm_depth} \
        -high_confidence_min_tumor_vaf {params.soft_vaf_cutoff} \
        -bqr_min_map_qual {params.min_map} -threads {threads} 2>&1 > {log}
        """

rule filter_sage:
    input:
        vcf = rules.run_sage.output.vcf
    output:
        vcf = config["cappseq_snv_pipeline"]["base_dir"] + "/02-vcfs/{sample}.sage.passed.vcf"
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        bcftools view -f PASS {input.vcf} -O vcf -o {output.vcf}
        """

# Flag positions with a high incidence of masked bases
rule flag_masked_pos:
    input:
        bam = lambda w: samplepaths[w.sample]
    output:
        bed_raw = temp(config["cappseq_snv_pipeline"]["base_dir"] + "/03-masked_pos/{sample}.maskedpos.bed"),
        bed = config["cappseq_snv_pipeline"]["base_dir"] + "/03-masked_pos/{sample}.maskedpos.bed.gz"
    params:
        script = os.path.join(config["cappseq_snv_pipeline"]["pipeline_dir"], "src/mask_n_sites.py"),
        n_threshold = config["cappseq_snv_pipeline"]["mask_threshold"],
        min_count = config["cappseq_snv_pipeline"]["mask_count"],
        panel_regions = config["cappseq_snv_pipeline"]["capture_space"]
    conda:
        "envs/bcftools.yaml"
    log:
        config["cappseq_snv_pipeline"]["base_dir"] + "/logs/{sample}.maskpos.log"
    shell:
        """
        {params.script} --input {input.bam} --regions {params.panel_regions} --output {output.bed_raw} --count {params.min_count} --fraction {params.n_threshold} > {log} &&
        bgzip -c {output.bed_raw} > {output.bed} && tabix -p bed {output.bed} >> {log}
        """

# Restrict to the captured regions, remove backlisted positions
rule restrict_to_capture:
    input:
        vcf = rules.filter_sage.output.vcf,
        bed = rules.flag_masked_pos.output.bed
    output:
        vcf = config["cappseq_snv_pipeline"]["base_dir"] + "/04-capturespace/{sample}.capspace.vcf"
    params:
        panel_regions = config["cappseq_snv_pipeline"]["capture_space"]
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        bedtools intersect -a {input.vcf} -header -b {params.panel_regions} | bedtools intersect -a - -header -b {input.bed} -v | awk -F '\\t' '$4 !~ /N/ && $5 !~ /N/' | bcftools norm -m +any -O vcf -o {output.vcf}
        """

# Generate review BAM files containing only the reads supporting these variants
rule review_consensus_reads:
    input:
        bam_cons = lambda w: samplepaths[w.sample],
        bam_uncons = lambda w: sample_uncons[w.sample],
        vcf = rules.restrict_to_capture.output.vcf
    output:
        tmp_sort = temp(config["cappseq_snv_pipeline"]["base_dir"] + "/06-supportingreads/{sample}/{sample}.umigrouped.sort.bam"),
        tmp_index = temp(config["cappseq_snv_pipeline"]["base_dir"] + "/06-supportingreads/{sample}/{sample}.umigrouped.sort.bam.bai"),
        consensus_bam = config["cappseq_snv_pipeline"]["base_dir"] + "/06-supportingreads/{sample}/{sample}.consensus.bam",
        grouped_bam = config["cappseq_snv_pipeline"]["base_dir"] + "/06-supportingreads/{sample}/{sample}.grouped.bam"
    threads: 2
    params:
        #fgbio = config["cappseq_snv_pipeline"]["fgbio_jar"],
        ref_genome = config["cappseq_snv_pipeline"]["ref_genome"],
        outdir = config["cappseq_snv_pipeline"]["base_dir"] + "/06-supportingreads/{sample}/"
    conda:
        "envs/fgbio.yaml"
    log:
        config["cappseq_snv_pipeline"]["base_dir"] + "/logs/{sample}.reviewconsensusvariant.log"
    shell:
        """
        samtools sort -@ 2 {input.bam_uncons} > {output.tmp_sort} && samtools index -@ 2 {output.tmp_sort} &&
        fgbio ReviewConsensusVariants --input {input.vcf} --consensus {input.bam_cons} --grouped-bam {output.tmp_sort} --ref {params.ref_genome} --output {params.outdir}/{wildcards.sample} --sample {wildcards.sample} 2>&1 > {log}
        """


rule vcf2maf_annotate:
    input:
        vcf = rules.restrict_to_capture.output.vcf
    output:
        vep_vcf = temp(config["cappseq_snv_pipeline"]["base_dir"] + "/04-capturespace/{sample}.capspace.vep.vcf"),
        maf = config["cappseq_snv_pipeline"]["base_dir"] + "/05-MAFs/{sample}.sage.maf"
    params:
        custom_enst = os.path.join(config["cappseq_snv_pipeline"]["pipeline_dir"], "resources/custom_enst.hg38.txt"),
        vep_data = config["cappseq_snv_pipeline"]["vep_data"],
        centre = config["cappseq_snv_pipeline"]["centre"],
        normal_name = lambda w: normal_ids[w.sample],
        ref_fasta = config["cappseq_snv_pipeline"]["ref_genome"],
        ref_ver = config["cappseq_snv_pipeline"]["ref_genome_ver"]
    threads: 2
    conda:
        "envs/vcf2maf.yaml"
    log:
        config["cappseq_snv_pipeline"]["base_dir"] + "/logs/{sample}.vcf2maf.log"
    shell:
        """
        vcf2maf.pl --input-vcf {input.vcf} --output-maf {output.maf} \
        --tumor-id {wildcards.sample} --normal-id {params.normal_name} \
        --vep-path $CONDA_PREFIX/bin/ \
        --vep-data {params.vep_data} --vep-forks {threads} \
        --custom-enst {params.custom_enst} --ref-fasta {params.ref_fasta} \
        --species homo_sapiens --ncbi-build {params.ref_ver} \
        --maf-center {params.centre} 2> {log} >> {log}
        """


rule augment_ssm:
    input:
        maf = rules.vcf2maf_annotate.output.maf,
        bam = lambda w: samplepaths[w.sample]
    output:
        maf = config["cappseq_snv_pipeline"]["base_dir"] + "/08-augmentssm/{sample}.sage.augment.maf"
    threads: 2
    conda:
        "envs/augment_ssm.yaml"
    params:
        script = os.path.join(config["cappseq_snv_pipeline"]["pipeline_dir"], "src/augment_ssm.py")
    shell: """
    python {params.script} \
    --bam {input.bam} \
    --maf {input.maf} \
    --output {output.maf} \
    --threads {threads}
    """


rule filter_alt_supp:
    input:
        maf = rules.augment_ssm.output.maf
    output:
        maf = config["cappseq_snv_pipeline"]["base_dir"] + "/09-altsupp/{sample}.sage.filtered.maf"
    params:
        min_alt_depth = config["cappseq_snv_pipeline"]["min_alt_depth"],
        min_germline_depth = config["cappseq_snv_pipeline"]["min_germline_depth"],
        tumor_af_diff = config["cappseq_snv_pipeline"]["tumor_af_diff"]
    shell: """
    awk -F '\t' '$42>={params.min_alt_depth} && $43>{params.min_germline_depth}' {input.maf} | \
    awk -F '\t' '$1 != "DNMT3A" || $42/$40 >= {params.tumor_af_diff} * $45/$43' | \
    grep -v ^HLA | grep -v ^FCGR > {output.maf}
    """

rule filter_repetitive_seq:
    """
    Remove mutations which are adjacent to a repetitive sequence.
    """
    input:
        maf = rules.filter_alt_supp.output.maf
    output:
        maf = config["cappseq_snv_pipeline"]["base_dir"] + "/10-filter_repeat/{sample}.sage.repeat_filt.maf"
    params:
        max_repeat_len = 6,
        ref_fasta = config["cappseq_snv_pipeline"]["ref_genome"]
    run:

        refseq = pyfaidx.Fasta(params.ref_fasta)
        with open(input.maf) as f, open(output.maf, "w") as o:
            for line in f:
                if line.startswith("#") or line.startswith("Hugo_Symbol"):
                    # Skip header lines.
                    o.write(line)
                    continue
                cols = line.split("\t")
                # Get alternate allele
                alt = cols[12]
                if alt == "-":
                    alt = cols[10]
                if len(alt) > 1:  # Handle large indels
                    alt = alt[0]
                # Get position of variant.
                
                pos = int(cols[6]) - 1  # Offset by 1 for MAF vs pyfaidx sequence.
                chrom = cols[4]

                # Determine length of repeat via reference sequence.
                repeat_len = 0
                while True:
                    pos += 1
                    base = refseq[chrom][pos]
                    if base != alt or repeat_len > params.max_repeat_len:
                        break
                    repeat_len += 1

                if repeat_len < params.max_repeat_len:
                    # Not a repeat (or not sufficiently long)
                    o.write(line)


# Filter output via GenomAD frequencies and position blacklist
rule filter_maf:
    input:
        maf = rules.filter_repetitive_seq.output.maf
    output:
        maf = config["cappseq_snv_pipeline"]["base_dir"] + "/99-final/{sample}.sage.blacklist.maf"
    params:
        exac_freq = float(config["cappseq_snv_pipeline"]["exac_max_freq"]),
        blacklist = os.path.join(config["cappseq_snv_pipeline"]["pipeline_dir"], "resources/capture-hg38.clean_blacklist.txt"),
        filter_cols = ["gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_ASJ_AF", "gnomAD_EAS_AF", "gnomAD_FIN_AF", "gnomAD_NFE_AF", "gnomAD_OTH_AF", "gnomAD_SAS_AF"]
    run:
        # Load blacklist
        blacklist_pos = []
        with open(params.blacklist) as f:
            for line in f:
                if line.startswith("#"):
                    continue  # Ignore comment lines
                # Assuming the first column of this file is chrom:pos
                line = line.rstrip("\n").rstrip("r")
                pos = line.split("\t")[0]
                blacklist_pos.append(pos)
        blacklist_pos = set(blacklist_pos)

        # Load variants
        in_maf = pd.read_csv(input.maf, sep ="\t", comment = "#")
        if in_maf.shape[0] == 0:  # Empty input MAF file
            in_maf.to_csv(output.maf, sep="\t", header=True, index = False)
        else:
            # Filter based on allele frequencies
            # This is slow and inefficient, but is probably fine for this implementation
            filtered_maf = in_maf
            for f_col in params.filter_cols:
                filtered_maf = filtered_maf[(numpy.isnan(filtered_maf[f_col])) | (filtered_maf[f_col] < params.exac_freq)]

            # Filter out any variants at positions overlapping the blacklist
            blacklist_maf = filtered_maf
            filtered_maf["key"] = filtered_maf["Chromosome"].astype(str) + ":" + filtered_maf["Start_Position"].astype(str)
            blacklist_maf = blacklist_maf[filtered_maf["key"].isin(blacklist_pos) == False]

            blacklist_maf.to_csv(output.maf, sep="\t", header=True, index = False)


# Generate IGV screenshots of these variants
rule igv_screenshot_variants:
    input:
        tbam = lambda w: samplepaths[w.sample],
        nbam = lambda w: normals[w.sample],
        maf = rules.filter_maf.output.maf
    output:
        html = config["cappseq_snv_pipeline"]["base_dir"] + "/07-IGV/{sample}_report.html"
    params:
        refgenome = config["cappseq_snv_pipeline"]["ref_genome"],
        cytoband = os.path.join(config['cappseq_snv_pipeline']['pipeline_dir'], config["cappseq_snv_pipeline"]["cytoband"]),
        genes =  os.path.join(config["cappseq_snv_pipeline"]['pipeline_dir'] , config["cappseq_snv_pipeline"]["genetrack"]),
        sample_name = lambda w: w.sample
    conda:
        "envs/igv.yaml"
    log:
        config["cappseq_snv_pipeline"]["base_dir"] + "/logs/{sample}.igv.log"
    shell:
        """
        create_report --fasta {params.refgenome} --type mutation --tracks {input.tbam} {input.nbam} {params.genes} --flanking 1500 --output {output.html} --standalone --title {params.sample_name} --ideogram {params.cytoband} {input.maf} > {log}
        """

rule all:
    input:
        expand(rules.filter_maf.output.maf, sample = samples),
        expand(rules.igv_screenshot_variants.output.html, sample=samples),
        expand(rules.review_consensus_reads.output.grouped_bam, sample=samples)
    default_target: True
