include: "tiling_common.smk"


rule tiling_extendRegions:
    conda:
        "../envs/default.yaml"
    input:
        regions=lambda wc: samples[wc.sample]["bed_file"],
        genome_file=config["reference"]["genome"],
    output:
        "results/tiling/{sample}/regions/extended.bed.gz",
    params:
        variant_edge_exclusion=config["tiling"]["variant_edge_exclusion"],
    log:
        "results/logs/tiling/extendRegions.{sample}.log",
    shell:
        """
        bedtools slop -i {input.regions} -g {input.genome_file} -b {params.variant_edge_exclusion} | \
        bgzip -c > {output} 2> {log}
        """


rule tiling_splitUpForStrategies:
    conda:
        "../envs/default.yaml"
    input:
        regions=lambda wc: tiling_getInputFile(wc.sample),
    output:
        noTiling="results/tiling/{sample}/strategies/regions.noTiling.bed.gz",
        noTiling_tmp=temp("results/tiling/{sample}/strategies/regions.noTiling.bed"),
        twoTiles="results/tiling/{sample}/strategies/regions.twoTiles.bed.gz",
        twoTiles_tmp=temp("results/tiling/{sample}/strategies/regions.twoTiles.bed"),
        centeredTiles="results/tiling/{sample}/strategies/regions.centeredTiles.bed.gz",
    params:
        oligo_length=config["oligo_length"],
        min_two_tiles=min(
            [
                config["tiling"]["strategies"]["two_tiles"]["min"],
                2 * config["oligo_length"] - 2 * config["tiling"]["min_overlap"],
            ]
        ),
    log:
        "results/logs/tiling/splitUpForStrategies.{sample}.log",
    shell:
        """
        zcat {input.regions} | \
        awk -v OFS="\\t" -v oligo_length={params.oligo_length} -v min_two_tiles={params.min_two_tiles} \
        '{{ 
            if ($3 - $2 <= oligo_length) {{ print $0 > {output.noTiling_tmp} }}
            else if ($3 - $2 <= min_two_tiles) {{ print $0 > {output.twoTiles_tmp} }}
            else {{ print $0 }}
        }}' | bgzip -c > {output.centeredTiles_tmp} 2> {log};
        cat {output.noTiling_tmp} | bgzip -c > {output.noTiling} 2>> {log};
        cat {output.twoTiles_tmp} | bgzip -c > {output.twoTiles} 2>> {log};
        """

rule tiling_strategy_twoTiles:
    """Tiling strategy: two tiles per region."""
    conda:
        "../envs/default.yaml"
    input:
        regions="results/tiling/{sample}/strategies/regions.twoTiles.bed.gz",
        genome_file=config["reference"]["genome"],
    output:
        "results/tiling/{sample}/strategies/regions.twoTiles.processed.bed.gz",
    params:
        oligo_length=config["oligo_length"],
        extension=config["tiling"]["variant_edge_exclusion"] if config["tiling"]["strategies"]["two_tiles"]["include_variant_edge"] else 0,
    shell:
        """
        bedtools slop -i {input.regions} -g {input.genome_file} -b {params.extension} | \
        awk -v OFS="\\t" -v oligo_length={params.oligo_length} \
        '{{ 
            start = $2;end = $3; $3 = start + oligo_length; print $0;
            $2 = end - oligo_length; $3 = end; print $0 
        }}' | \
        bgzip -c > {output} 2> {log}
        """

rule tiling_strategy_centeredTiles:
    """Tiling strategy: two tiles per region."""
    conda:
        "../envs/default.yaml"
    input:
        regions="results/tiling/{sample}/strategies/regions.centeredTiles.bed.gz",
        genome_file=config["reference"]["genome"],
    output:
        "results/tiling/{sample}/strategies/regions.centeredTiles.processed.bed.gz",
    params:
        oligo_length=config["oligo_length"],
        extension=config["tiling"]["variant_edge_exclusion"] if config["tiling"]["strategies"]["two_tiles"]["include_variant_edge"] else 0,
    shell:
        """
        bedtools slop -i {input.regions} -g {input.genome_file} -b {params.extension} | \
        awk -v OFS="\\t" -v oligo_length={params.oligo_length} \
        '{{ 
            start = $2;end = $3; $3 = start + oligo_length; print $0;
            $2 = end - oligo_length; $3 = end; print $0 
        }}' | \
        bgzip -c > {output} 2> {log}
        """
