
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
