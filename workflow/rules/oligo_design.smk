
rule oligo_design_getSequencesInclVariants:
    """Rule to map variants to designed regions, 
    filter not matched variants, 
    removes not used regions (optinal), 
    and generate ref/alt sequences"""
    conda:
        "../envs/oligo_design.yaml"
    input:
        regions="results/tiling/{sample}/regions.tiles.bed.gz",
        variants=lambda wc: samples.loc[wc.sample]["vcf_file"],
        ref=config["reference"]["fasta"],
        genome_file=config["reference"]["genome"],
        script=getScript("oligo_design/getSequencesInclVariants.py"),
    output:
        variants="results/oligo_design/{sample}/design.variants.vcf.gz",
        variants_removed="results/oligo_design/{sample}/removed.variants.vcf.gz",
        regions="results/oligo_design/{sample}/design.regions.bed.gz",
        regions_removed="results/oligo_design/{sample}/removed.regions.bed.gz",
        design="results/oligo_design/{sample}/design.fa",
        design_map="results/oligo_design/{sample}/design.variant_region_map_map.tsv.gz",
    params:
        variant_edge_exclusion=config["tiling"]["variant_edge_exclusion"],
        use_most_centered_region="--use-most-centered-region-for-variant"
        if config["oligo_design"]["variants"]["use_most_centered_region"]
        else "--use-all-regions-for-variant",
        remove_regions_without_variants="--remove-regions-without-variants"
        if config["oligo_design"]["variants"]["remove_unused_regions"]
        else "--keep-regions-without-variants",
    log:
        "logs/oligo_design/getSequencesInclVariants.{sample}.log",
    shell:
        """
        python {input.script} \
        --input-regions {input.regions} \
        --input-variants {input.variants} \
        --output-variants {output.variants} \
        --output-variants-removed {output.variants_removed} \
        --output-regions {output.regions} \
        --output-regions-without-variants {output.regions_removed} \
        --output-design {output.design} \
        --output-design-map {output.design_map} \
        --reference {input.ref} \
        {params.remove_regions_without_variants} --variant-edge-exclusion {params.variant_edge_exclusion} {params.use_most_centered_region}
        """

rule oligo_design_filterOligos:
    """Remove oligos overlapping with TSS, CTCF, 
    restriction sites and simple repeats, 
    as well as oligos with too many homopolymers. 
    """
    conda:
        "../envs/filter.yaml"
    input:
        regions="results/oligo_design/{sample}/design.regions.bed.gz",
        design="results/oligo_design/{sample}/design.fa",
        script=getScript("oligo_design/filterOligos.py"),
        path="results/oligo_design/{sample}/"
    output:
        regions="results/oligo_design/{sample}/filtered.regions.bed.gz",
        design="results/oligo_design/{sample}/filtered.design.fa",
        variant_ids="results/oligo_design/{sample}/filtered.var.ids.txt",
    params:
        repeats=config["oligo_design"]["filtering"]["max_simple_repeat_fraction"], 
        max_hom=config["oligo_design"]["filtering"]["max_homopolymer_length"],
    shell:
        """
        python {input.script} \
        --seqs {input.design} \
        --repeat {params.repeats} \
        --max_homopolymer_length {params.max_hom} \
        --outpath {input.path}
        """

rule filter_variants:
    """Retain only those variants which are still included after filtering the oligos.
    """
    input:
        filtered_ids="results/oligo_design/{sample}/filtered.var.ids.txt",
        variants="results/oligo_design/{sample}/design.variants.vcf.gz"
    output:
        filtered_variants="results/oligo_design/{sample}/filtered.variants.vcf.gz"
    shell:
        """
         cat <(zgrep '^#' {input.variants} ) <( awk 'NR=FNR{{a[$3][$0]}} $0 in a {{for (i in a[$0]) print i}}' \
         <(zcat {input.variants}) {input.filtered_ids}) | bgzip -c > {output.filtered_variants}
        #"""
    