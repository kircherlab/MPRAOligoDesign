##################
#### variants ####
##################


rule oligo_design_getSequencesInclVariants:
    """
    Rule to map variants to designed regions, 
    filter not matched variants, 
    removes not used regions (optional), 
    and generate ref/alt sequences
    """
    conda:
        "../envs/oligo_design.yaml"
    input:
        regions=lambda wc: "results/tiling/{sample}/regions.tiles.bed.gz"
        if isVariantsAndRegionsSample(wc.sample)
        else "results/centering/{sample}/regions.centered.bed.gz",
        variants=lambda wc: datasets.loc[wc.sample]["vcf_file"],
        ref=config["reference"]["fasta"],
        genome_file=config["reference"]["genome"],
        script=getScript("oligo_design/getSequencesInclVariants.py"),
    output:
        variants="results/oligo_design/{sample}/design_variants.variants.vcf.gz",
        variants_removed="results/oligo_design/{sample}/removed.variants.vcf.gz",
        regions="results/oligo_design/{sample}/design_variants.regions.bed.gz",
        regions_removed="results/oligo_design/{sample}/removed.regions.bed.gz",
        design="results/oligo_design/{sample}/design_variants.fa",
        design_variant_map="results/oligo_design/{sample}/design_variants.variant_region_map.tsv.gz",
        design_region_map="results/oligo_design/{sample}/design_variants.region_map.tsv.gz",
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
        --output-design-variant-map {output.design_variant_map} \
        --output-design-region-map {output.design_region_map} \
        --reference {input.ref} \
        {params.remove_regions_without_variants} --variant-edge-exclusion {params.variant_edge_exclusion} {params.use_most_centered_region} &> {log}
        """


rule oligo_design_variants_filterOligos:
    """
    Remove oligos overlapping with TSS, CTCF, 
    restriction sites and simple repeats, 
    as well as oligos with too many homopolymers. 
    """
    conda:
        "../envs/filter.yaml"
    input:
        regions="results/oligo_design/{sample}/design_variants.regions.bed.gz",
        design="results/oligo_design/{sample}/design_variants.fa",
        variant_map="results/oligo_design/{sample}/design_variants.variant_region_map.tsv.gz",
        region_map="results/oligo_design/{sample}/design_variants.region_map.tsv.gz",
        simple_repeats=getReference("simpleRepeat.bed.gz"),
        tss=getReference("TSS_pos.bed.gz"),
        ctcf=getReference("CTCF-MA0139-1_intCTCF_fp25.hg38.bed.gz"),
        script=getScript("oligo_design/filterOligos.py"),
    output:
        out_variant_map="results/oligo_design/{sample}/design_variants_filtered.variant_region_map.tsv.gz",
        out_region_map="results/oligo_design/{sample}/design_variants_filtered.region_map.tsv.gz",
        statistic="results/oligo_design/{sample}/design_variants_filter.log",
    params:
        repeats=config["oligo_design"]["filtering"]["max_simple_repeat_fraction"],
        max_hom=config["oligo_design"]["filtering"]["max_homopolymer_length"],
        remove_regions_without_variants="--remove-regions-without-variants"
        if config["oligo_design"]["variants"]["remove_unused_regions"]
        else "--keep-regions-without-variants",
    log:
        "logs/oligo_design/filterOligos.{sample}.log",
    shell:
        """
        python {input.script} \
        --seqs {input.design} \
        --regions {input.regions} \
        --variant-map {input.variant_map} \
        --map {input.region_map} \
        --repeat {params.repeats} \
        --max_homopolymer_length {params.max_hom} \
        --tss-positions {input.tss} \
        --ctcf-motifs {input.ctcf} \
        --simple-repeats {input.simple_repeats} \
        {params.remove_regions_without_variants} \
        --output-map {output.out_region_map} \
        --output-variant-map {output.out_variant_map} \
        > {output.statistic} 2> {log}
        """


rule oligo_design_variants_filter_regions:
    """Filter the bed file using the filtered map."""
    conda:
        "../envs/default.yaml"
    input:
        map="results/oligo_design/{sample}/design_variants_filtered.region_map.tsv.gz",
        regions="results/oligo_design/{sample}/design_variants.regions.bed.gz",
    output:
        regions="results/oligo_design/{sample}/design_variants_filtered.regions.bed.gz",
    log:
        "logs/oligo_design/filter_regions.{sample}.log",
    shell:
        """
        awk 'NR=FNR{{a[$4][$0]}} $1 in a {{for (i in a[$1]) print i}}' \
        <(zcat {input.regions}) <(zcat {input.map}) | sort -k1,1 -k2,2n | uniq | \
        bgzip -c > {output.regions} 2> {log}"""


rule oligo_design_variants_filter_seqs:
    """Filter the bed file using the filtered map."""
    conda:
        "../envs/default.yaml"
    input:
        map="results/oligo_design/{sample}/design_variants_filtered.region_map.tsv.gz",
        seqs="results/oligo_design/{sample}/design_variants.fa",
    output:
        seqs="results/oligo_design/{sample}/design_variants_filtered.design.fa",
    log:
        "logs/oligo_design/filter_seqs.{sample}.log",
    shell:
        """
        awk 'NR=FNR {{a[$2]}} {{if ($1 ~ /^>/) id=substr($1,2); if (id in a) print $0}}' \
        <(zcat {input.map}) {input.seqs} > {output.seqs} 2> {log}
        """


rule oligo_design_variants_filter_variants:
    """
    Retain only those variants which are still included after filtering the oligos.
    """
    conda:
        "../envs/default.yaml"
    input:
        map="results/oligo_design/{sample}/design_variants_filtered.variant_region_map.tsv.gz",
        variants="results/oligo_design/{sample}/design_variants.variants.vcf.gz",
    output:
        filtered_variants="results/oligo_design/{sample}/design_variants_filtered.variants.vcf.gz",
    log:
        "logs/oligo_design/filter_variants.{sample}.log",
    shell:
        """
        (
            zgrep '^#' {input.variants};
            awk 'NR=FNR{{a[$3][$0]}} $1 in a {{for (i in a[$1]) print i}}' \
            <(zcat {input.variants}) <(zcat {input.map}) | sort -k1,1 -k2,2n -k3,3 | uniq;
        ) | \
        bgzip -c > {output.filtered_variants} 2> {log}
        """


#####################
#### Region only ####
#####################
rule oligo_design_regions_getSequences:
    """
    retrieve the sequences for the design using only regions
    """
    conda:
        "../envs/default.yaml"
    input:
        regions=lambda wc: datasets.loc[wc.sample, "bed_file"],
        ref=config["reference"]["fasta"],
    output:
        design="results/oligo_design/{sample}/design_regions.fa",
    params:
        region_size=config["oligo_length"],
    log:
        "logs/oligo_design/regions_getSequences.{sample}.log",
    shell:
        """
        error=`zcat {input.regions} | awk '{{if ($3-$2 != {params.region_size}) {{print $4}}}}'`
        if [ -n "${{error}}" ]; then
            echo "ERROR: The following regions have a different size than the oligo length: $error"
            exit 1
        fi
        bedtools getfasta -s -nameOnly -fi {input.ref} -bed {input.regions} | \
        sed 's/(+)$//' | sed 's/(-)$//' > {output} 2> {log}
        """


rule oligo_design_regions_getRegionMap:
    """
    retrieve the sequences for the design using only regions
    """
    conda:
        "../envs/default.yaml"
    input:
        regions=lambda wc: datasets.loc[wc.sample, "bed_file"],
        design="results/oligo_design/{sample}/design_regions.fa",
    output:
        design_map="results/oligo_design/{sample}/design_regions.region_map.tsv.gz",
    log:
        "logs/oligo_design/regions_getRegionMap.{sample}.log",
    shell:
        """
        (
            echo -e "Region\\tID";
            paste <(zcat {input.regions} | cut -f 4) <(cat {input.design} | egrep "^>" | sed 's/^>//')
        ) | gzip -c > {output} 2> {log} 
        """


rule oligo_design_regions_filterOligos:
    """
    Remove oligos overlapping with TSS, CTCF, 
    restriction sites and simple repeats, 
    as well as oligos with too many homopolymers. 
    """
    conda:
        "../envs/filter.yaml"
    input:
        regions=lambda wc: datasets.loc[wc.sample, "bed_file"],
        design="results/oligo_design/{sample}/design_regions.fa",
        design_map="results/oligo_design/{sample}/design_regions.region_map.tsv.gz",
        simple_repeats=getReference("simpleRepeat.bed.gz"),
        tss=getReference("TSS_pos.bed.gz"),
        ctcf=getReference("CTCF-MA0139-1_intCTCF_fp25.hg38.bed.gz"),
        script=getScript("oligo_design/filterOligos.py"),
    output:
        out_map="results/oligo_design/{sample}/design_regions_filtered.region_map.tsv.gz",
        statistic="results/oligo_design/{sample}/design_regions_filter.log",
    params:
        repeats=config["oligo_design"]["filtering"]["max_simple_repeat_fraction"],
        max_hom=config["oligo_design"]["filtering"]["max_homopolymer_length"],
    log:
        "logs/oligo_design/regions_filterOligos.{sample}.log",
    shell:
        """
        python {input.script} \
        --seqs {input.design} \
        --regions {input.regions} \
        --map {input.design_map} \
        --repeat {params.repeats} \
        --max_homopolymer_length {params.max_hom} \
        --tss-positions {input.tss} \
        --ctcf-motifs {input.ctcf} \
        --simple-repeats {input.simple_repeats} \
        --output-map {output.out_map} \
        > {output.statistic} 2> {log}
        """


rule oligo_design_regions_filter_regions:
    """Filter the bed file using the filtered map."""
    conda:
        "../envs/default.yaml"
    input:
        map="results/oligo_design/{sample}/design_regions_filtered.region_map.tsv.gz",
        regions=lambda wc: datasets.loc[wc.sample, "bed_file"],
    output:
        regions="results/oligo_design/{sample}/design_regions_filtered.regions.bed.gz",
    log:
        "logs/oligo_design/regions_filter_regions.{sample}.log",
    shell:
        """
        awk 'NR=FNR{{a[$4][$0]}} $2 in a {{for (i in a[$1]) print i}}' \
        <(zcat {input.regions}) <(zcat {input.map}) | sort -k1,1 -k2,2n | uniq | \
        bgzip -c > {output.regions} 2> {log}
        """


rule oligo_design_regions_filter_seqs:
    """Filter the bed file using the filtered map."""
    conda:
        "../envs/default.yaml"
    input:
        map="results/oligo_design/{sample}/design_regions_filtered.region_map.tsv.gz",
        seqs="results/oligo_design/{sample}/design_regions.fa",
    output:
        seqs="results/oligo_design/{sample}/design_regions_filtered.design.fa",
    log:
        "logs/oligo_design/filter_seqs.{sample}.log",
    shell:
        """
        awk 'NR=FNR {{a[$2]}} {{if ($1 ~ /^>/) id=substr($1,2); if (id in a) print $0}}' \
        <(zcat {input.map}) {input.seqs} > {output.seqs} 2> {log}
        """


#####################
#### Sequence only ####
#####################


rule oligo_design_seq_getsequenceMap:
    """
    retrieve the sequences for the design using only regions
    """
    conda:
        "../envs/default.yaml"
    input:
        design=lambda wc: datasets.loc[wc.sample, "fasta_file"],
    output:
        design_map="results/oligo_design/{sample}/design_sequences.sequence_map.tsv.gz",
    params:
        region_size=config["oligo_length"],
    log:
        "logs/oligo_design/seq_getsequenceMap.{sample}.log",
    shell:
        """
        error=`cat {input.design} | sed 's/\\r//' | awk -v "OFS=\\t" '{{if($0 ~ /^>/ ) {{id=$0}} else if(length($1) != {params.region_size}) {{print id, length($1)}}}}'`
        if [ -n "${{error}}" ]; then
            echo "ERROR: The following sequence have a different size than the oligo length: ${{error}}"
            exit 1
        fi
        (
            echo -e "ID";
            cat {input.design} | sed 's/\\r//' | egrep "^>" | sed 's/^>//';
        ) | \
        gzip -c > {output} 2> {log} 
        """


rule oligo_design_seq_filterOligos:
    """
    Remove oligos overlapping with TSS, CTCF, 
    restriction sites and simple repeats, 
    as well as oligos with too many homopolymers. 
    """
    conda:
        "../envs/filter.yaml"
    input:
        design=lambda wc: datasets.loc[wc.sample, "fasta_file"],
        design_map="results/oligo_design/{sample}/design_sequences.sequence_map.tsv.gz",
        simple_repeats=getReference("simpleRepeat.bed.gz"),
        tss=getReference("TSS_pos.bed.gz"),
        ctcf=getReference("CTCF-MA0139-1_intCTCF_fp25.hg38.bed.gz"),
        script=getScript("oligo_design/filterOligos.py"),
    output:
        out_map="results/oligo_design/{sample}/design_sequences_filtered.sequence_map.tsv.gz",
        statistic="results/oligo_design/{sample}/design_sequences_filter.log",
    params:
        repeats=config["oligo_design"]["filtering"]["max_simple_repeat_fraction"],
        max_hom=config["oligo_design"]["filtering"]["max_homopolymer_length"],
    log:
        "logs/oligo_design/seq_filterOligos.{sample}.log",
    shell:
        """
        python {input.script} \
        --seqs {input.design} \
        --map {input.design_map} \
        --repeat {params.repeats} \
        --max_homopolymer_length {params.max_hom} \
        --tss-positions {input.tss} \
        --ctcf-motifs {input.ctcf} \
        --simple-repeats {input.simple_repeats} \
        --output-map {output.out_map} \
        > {output.statistic} 2> {log}
        """


rule oligo_design_sequences_filter_seqs:
    """Filter the bed file using the filtered map."""
    conda:
        "../envs/default.yaml"
    input:
        map="results/oligo_design/{sample}/design_sequences_filtered.sequence_map.tsv.gz",
        seqs=lambda wc: datasets.loc[wc.sample, "fasta_file"],
    output:
        seqs="results/oligo_design/{sample}/design_sequences_filtered.design.fa",
    log:
        "logs/oligo_design/filter_seqs.{sample}.log",
    shell:
        """
        awk 'NR=FNR {{a[$1]}} {{if ($1 ~ /^>/) id=substr($1,2); if (id in a) print $0}}' \
        <(zcat {input.map}) {input.seqs} > {output.seqs} 2> {log}
        """
