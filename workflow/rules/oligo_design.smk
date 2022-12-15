
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
        design_map="results/oligo_design/{sample}/design.variant_region_map.tsv.gz",
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
        {params.remove_regions_without_variants} --variant-edge-exclusion {params.variant_edge_exclusion} {params.use_most_centered_region} &> {log}
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
        seq_map="results/oligo_design/{sample}/design.variant_region_map.tsv.gz",
        simple_repeats=getReference("simpleRepeat.bed.gz"),
        tss=getReference("TSS_pos.bed.gz"),
        ctcf=getReference("CTCF-MA0139-1_intCTCF_fp25.hg38.bed.gz"),
        script=getScript("oligo_design/filterOligos.py"),
    output:
        #design="results/oligo_design/{sample}/filtered.design.fa",
        out_map="results/oligo_design/{sample}/filtered.variant_region_map.tsv.gz",
        #variant_ids="results/oligo_design/{sample}/filtered.var.ids.txt",
        vatiants_failed="results/oligo_design/{sample}/failed.var.ids.txt",
        statistic="results/oligo_design/{sample}/filter.log",
    params:
        repeats=config["oligo_design"]["filtering"]["max_simple_repeat_fraction"],
        max_hom=config["oligo_design"]["filtering"]["max_homopolymer_length"],
    log:
        "logs/oligo_design/filterOligos.{sample}.log",
    shell:
        """
        python {input.script} \
        --seqs {input.design} \
        --regions {input.regions} \
        --seq-map {input.seq_map} \
        --repeat {params.repeats} \
        --max_homopolymer_length {params.max_hom} \
        --tss-positions {input.tss} \
        --ctcf-motifs {input.ctcf} \
        --simple-repeats {input.simple_repeats} \
        --output-map {output.out_map} \
        --output-failed-variants {output.vatiants_failed} \
        > {output.statistic} 2> {log}
        """


rule oligo_design_filter_regions:
    """Filter the bed file using the filtered map."""
    conda:
        "../envs/default.yaml"
    input:
        map="results/oligo_design/{sample}/filtered.variant_region_map.tsv.gz",
        regions="results/oligo_design/{sample}/design.regions.bed.gz",
    output:
        regions="results/oligo_design/{sample}/filtered.regions.bed.gz",
    log:
        "logs/oligo_design/filter_regions.{sample}.log",
    shell:
        """
        awk 'NR=FNR{{a[$4][$0]}} $2 in a {{for (i in a[$2]) print i}}' \
        <(zcat {input.regions}) <(zcat {input.map}) | \
        bgzip -c > {output.regions} 2> {log}"""


rule oligo_design_filter_seqs:
    """Filter the bed file using the filtered map."""
    conda:
        "../envs/default.yaml"
    input:
        map="results/oligo_design/{sample}/filtered.variant_region_map.tsv.gz",
        seqs="results/oligo_design/{sample}/design.fa",
    output:
        seqs="results/oligo_design/{sample}/filtered.design.fa",
    log:
        "logs/oligo_design/filter_seqs.{sample}.log",
    shell:
        """
        awk 'NR=FNR {{a[$3]; a[$4]}} {{if ($1 ~ /^>/) id=substr($1,2); if (id in a) print $0}}' \
        <(zcat {input.map}) {input.seqs} > {output.seqs} 2> {log}"""


rule oligo_design_filter_variants:
    """Retain only those variants which are still included after filtering the oligos.
    """
    conda:
        "../envs/default.yaml"
    input:
        map="results/oligo_design/{sample}/filtered.variant_region_map.tsv.gz",
        variants="results/oligo_design/{sample}/design.variants.vcf.gz",
    output:
        filtered_variants="results/oligo_design/{sample}/filtered.variants.vcf.gz",
    log:
        "logs/oligo_design/filter_variants.{sample}.log",
    shell:
        """
        cat <(zgrep '^#' {input.variants} ) <( awk 'NR=FNR{{a[$3][$0]}} $1 in a {{for (i in a[$1]) print i}}' \
        <(zcat {input.variants}) <(zcat {input.map})) | \
        bgzip -c > {output.filtered_variants} 2> {log}
        #"""


rule oligo_design_add_adapters:
    """Add adapters to the final designed sequences.
    """
    conda:
        "../envs/default.yaml"
    input:
        designs="results/oligo_design/{sample}/filtered.design.fa",
    output:
        designs="results/oligo_design/{sample}/final.design.fa",
    log:
        "logs/oligo_design/add_adapters.{sample}.log",
    params:
        left=config["oligo_design"]["adapters"]["left"],
        right=config["oligo_design"]["adapters"]["right"],
    shell:
        """
        sed -e '/^>/! s/^/{params.left}/' -e '/^>/!  s/$/{params.right}/' {input.designs} > \
        {output.designs} 2> {log}
        """
