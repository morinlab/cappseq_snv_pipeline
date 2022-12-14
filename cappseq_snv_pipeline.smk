#!/usr/bin/env snakemake

import os
import pandas
import numpy
import snakemake

snakemake.utils.min_version("7")

# Load config file
configfile: "config/cappseq_snv_pipeline.yaml"
configpath = "config/cappseq_snv_pipeline.yaml"

# Check that the config file has all the required parameters
pathkeys = {"samplelist", "basedir", "ref_genome", "hotspots_vcf", "capture_space", "ensembl", "unmatched_normal", "custom_enst", "vep_data", "fgbio_jar"}
for ckey, attribute in config["cappseq_snv_pipeline"].items():
    if attribute == "__UPDATE__":
        # Placeholder value in config. Warn user
        raise AttributeError(f"\'__UPDATE__\' found for \'{ckey}\' in config file \'{configpath}\'. Please ensure the config file is updated with parameters relevant for your analysis")
    # Check that required filepaths exist
    if ckey in pathkeys:
        if not os.path.exists(attribute):
            raise AttributeError(f"Unable to locate \'{attribute}\' for key \'{ckey}\' in config file \'{configpath}\': No such file or directory")


# Load input files
samplelist = config["cappseq_snv_pipeline"]["samplelist"]

samples = []
samplepaths = {}
sample_uncons = {}

with open(samplelist) as f:
    i = 0
    for line in f:
        # Assuming a two-column sample file, with sample name as column 1, and
        # a BAM/CRAM path as column 2
        i += 1
        # Ignore comment lines
        if line.startswith("#"):
            continue
        line = line.rstrip("\n").rstrip("\r")
        cols = line.split("\t")
        try:
            sample_name = cols[0]
            sample_path = cols[1]
            unmerged_path = cols[2]
        except IndexError as e:
            raise AttributeError(f"Unable to parse line {i} of sample file \'{samplelist}\': Expected three columns specifying sample name, consensus BAM path, and pre-consensus BAM path") from e
        # Sanity check that the filepath exists
        if not os.path.exists(sample_path):
            raise AttributeError(f"Unable to locate BAM/CRAM file \'{sample_path}\' for sample \'{sample_name}\': No such file or directory")
        # Store these samples
        if sample_name in samples:
            raise AttributeError(f"Duplicate sample name \'{sample_name}\' detected in sample file \'{samplelist}\'")
        samples.append(sample_name)
        samplepaths[sample_name] = sample_path
        sample_uncons[sample_name] = unmerged_path

# Run variant calling
rule run_sage:
    input:
        bam = lambda w: samplepaths[w.sample]
    output:
        vcf = config["cappseq_snv_pipeline"]["base_dir"] + "/01-SAGE/{sample}/{sample}.sage.vcf"
    params:
        ref_genome = config["cappseq_snv_pipeline"]["ref_genome"],
        ref_genome_version = "38" if config["cappseq_snv_pipeline"]["ref_genome_ver"] == "GRCh38" else "37",  # Should make more robust
        # Panel regions
        hotspots_vcf = config["cappseq_snv_pipeline"]["hotspots_vcf"],
        panel_regions = config["cappseq_snv_pipeline"]["capture_space"],
        ensembl = config["cappseq_snv_pipeline"]["ensembl"],
        # Normal
        normal_bam = config["cappseq_snv_pipeline"]["unmatched_normal"],
        normal_name = config["cappseq_snv_pipeline"]["normal_name"],
        # Miscellaneous
        max_depth = config["cappseq_snv_pipeline"]["max_depth"],
        min_map = config["cappseq_snv_pipeline"]["min_map_qual"]
    threads: 4
    log:
        config["cappseq_snv_pipeline"]["base_dir"] + "/logs/{sample}.sage_run.log"
    conda:
        "envs/sage.yaml"
    shell:
        """SAGE -tumor {wildcards.sample} -tumor_bam {input.bam} \
        -out {output.vcf} -ref_genome {params.ref_genome} \
        -ref_genome_version {params.ref_genome_version} \
        -hotspots {params.hotspots_vcf} -panel_bed {params.panel_regions} \
        -high_confidence_bed {params.panel_regions} \
        -ensembl_data_dir {params.ensembl} \
        -reference_bam {params.normal_bam} -reference {params.normal_name} \
        -max_read_depth {params.max_depth} -min_map_quality {params.min_map} \
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
        script = config["cappseq_snv_pipeline"]["mask_script"],
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
        bedtools intersect -a {input.vcf} -header -b {params.panel_regions} | bedtools intersect -a - -header -b {input.bed} -v | bcftools filter -e 'ALT=\"N\"' -O vcf | bcftools norm -m +any -O vcf -o {output.vcf}
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
        fgbio = config["cappseq_snv_pipeline"]["fgbio_jar"],
        ref_genome = config["cappseq_snv_pipeline"]["ref_genome"],
        outdir = config["cappseq_snv_pipeline"]["base_dir"] + "/06-supportingreads/{sample}/"
    conda:
        "envs/fgbio.yaml"
    log:
        config["cappseq_snv_pipeline"]["base_dir"] + "/logs/{sample}.reviewconsensusvariant.log"
    shell:
        """
        samtools sort -@ 2 {input.bam_uncons} > {output.tmp_sort} && samtools index -@ 2 {output.tmp_sort} &&
        java -jar {params.fgbio} ReviewConsensusVariants --input {input.vcf} --consensus {input.bam_cons} --grouped-bam {output.tmp_sort} --ref {params.ref_genome} --output {params.outdir}/{wildcards.sample} --sample {wildcards.sample} 2>&1 > {log}
        """


rule vcf2maf_annotate:
    input:
        vcf = rules.restrict_to_capture.output.vcf
    output:
        vep_vcf = temp(config["cappseq_snv_pipeline"]["base_dir"] + "/04-capturespace/{sample}.capspace.vep.vcf"),
        maf = config["cappseq_snv_pipeline"]["base_dir"] + "/05-MAFs/{sample}.sage.maf"
    params:
        custom_enst = config["cappseq_snv_pipeline"]["custom_enst"],
        vep_data = config["cappseq_snv_pipeline"]["vep_data"],
        centre = config["cappseq_snv_pipeline"]["centre"],
        normal_name = config["cappseq_snv_pipeline"]["normal_name"],
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

# Filter output via GenomAD frequencies and position blacklist
rule filter_maf:
    input:
        maf = rules.vcf2maf_annotate.output.maf
    output:
        maf = config["cappseq_snv_pipeline"]["base_dir"] + "/99-final/{sample}.sage.filtered.maf"
    params:
        exac_freq = float(config["cappseq_snv_pipeline"]["exac_max_freq"]),
        blacklist = config["cappseq_snv_pipeline"]["blacklist"],
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
        in_maf = pandas.read_csv(input.maf, sep ="\t", comment = "#")
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
        full_bam = lambda w: samplepaths[w.sample],
        consensus_bam = rules.review_consensus_reads.output.consensus_bam,
        grouped_bam = rules.review_consensus_reads.output.grouped_bam,
        maf = rules.filter_maf.output.maf
    output:
        outdir = directory(config["cappseq_snv_pipeline"]["base_dir"] + "/07-IGV/{sample}/")
    params:
        igv_script = config["cappseq_snv_pipeline"]["igv_script"],
        genome_build = "hg38" if config["cappseq_snv_pipeline"]["ref_genome_ver"] == "GRCh38" else "37"
    conda:
        "envs/igv.yaml"
    log:
        config["cappseq_snv_pipeline"]["base_dir"] + "/logs/{sample}.igv.log"
    shell:
        """
        python {params.igv_script} -m {input.maf} -b {input.full_bam} {input.consensus_bam} {input.grouped_bam} --genome_build {params.genome_build} -o {output.outdir} 2> {log}
        """

rule all:
    input:
        expand(rules.igv_screenshot_variants.output.outdir, sample = samples)
    default_target: True
