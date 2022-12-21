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
        noTiling="results/tiling/{sample}/strategies/regions.noTiles.bed.gz",
        noTiling_tmp=temp("results/tiling/{sample}/strategies/regions.noTiles.bed"),
        twoTiles="results/tiling/{sample}/strategies/regions.twoTiles.bed.gz",
        twoTiles_tmp=temp("results/tiling/{sample}/strategies/regions.twoTiles.bed"),
        centeredTiles="results/tiling/{sample}/strategies/regions.centeredTiles.bed.gz",
    params:
        oligo_length=config["oligo_length"],
        max_two_tiles=tiling_getMinTwoTiles(),
    log:
        "results/logs/tiling/splitUpForStrategies.{sample}.log",
    shell:
        """
        touch {output.noTiling_tmp};
        touch {output.twoTiles_tmp};
        (
            if [[ $(file -b --mime-type "{input.regions}") == 'application/gzip' ]]; then
                zcat {input.regions}
            else
                cat {input.regions}
            fi
        ) | \
        awk -v FS="\\t" -v OFS="\\t" \
        '{{ 
            if ($3 - $2 <= {params.oligo_length}) {{print > "{output.noTiling_tmp}" }}
            else if ($3 - $2 <= {params.max_two_tiles}) {{ print > "{output.twoTiles_tmp}" }}
            else {{ print }}
        }}' | bgzip -c > {output.centeredTiles} 2> {log};
        cat {output.noTiling_tmp} | bgzip -c > {output.noTiling} 2>> {log};
        cat {output.twoTiles_tmp} | bgzip -c > {output.twoTiles} 2>> {log};
        """


rule tiling_strategy_noTiles:
    """Tiling strategy: no tiles per region."""
    conda:
        "../envs/default.yaml"
    input:
        regions="results/tiling/{sample}/strategies/regions.noTiles.bed.gz",
    output:
        "results/tiling/{sample}/strategies/regions.noTiles.processed.bed.gz",
    params:
        oligo_length=config["oligo_length"],
    log:
        "results/logs/tiling/strategy_noTiles.{sample}.log",
    shell:
        """
        bgzip -dc {input.regions} | \
        awk -v OFS="\\t" -v oligo_length={params.oligo_length} \
        'function ceil(x){{return int(x)+(x>int(x))}};
        function floor(x, y){{y=int(x); return(x<y?y-1:y)}};
        {{ 
            strand="none";
            if ($6 == "+") {{strand="fwd"}}
            else if ($6 == "-") {{strand="rev"}};
            missing = (oligo_length - ($3 - $2))/2;
            $2 = $2-ceil(missing);
            $3 = $3+floor(missing);
            $4=$4"_"strand"_tile1-1"
            print;
        }}' | \
        bgzip -c > {output} 2> {log}
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
        extension=config["tiling"]["variant_edge_exclusion"]
        if config["tiling"]["strategies"]["two_tiles"]["include_variant_edge"]
        else 0,
    log:
        "results/logs/tiling/strategy_twoTiles.{sample}.log",
    shell:
        """
        bedtools slop -i {input.regions} -g {input.genome_file} -b {params.extension} | \
        awk -v FS="\\t" -v OFS="\\t" -v oligo_length={params.oligo_length} \
        '{{ 
            strand="none";
            if ($6 == "+") {{strand="fwd"}}
            else if ($6 == "-") {{strand="rev"}};
            start = $2; end = $3; id = $4;
            $3 = start + oligo_length; $4=id"_"strand"_tile1-2"; print $0;
            $2 = end - oligo_length; $3 = end; $4=id"_"strand"_tile2-2"; print $0;
        }}' | \
        bgzip -c > {output} 2> {log}
        """


rule tiling_strategy_centeredTiles:
    """Tiling strategy: two tiles per region."""
    conda:
        "../envs/default.yaml"
    input:
        regions="results/tiling/{sample}/strategies/regions.centeredTiles.bed.gz",
        script=getScript("tiling/centerTiling.py"),
    output:
        "results/tiling/{sample}/strategies/regions.centeredTiles.processed.bed.gz",
    params:
        oligo_length=config["oligo_length"],
        min_overlap=config["tiling"]["min_overlap"],
    log:
        "results/logs/tiling/strategy_centeredTiles.{sample}.log",
    shell:
        """
        python {input.script} \
        --input {input.regions} \
        --oligo-length {params.oligo_length} \
        --min-overlap {params.min_overlap} \
        --output >( bgzip -c > {output}) &> {log}
        """


rule tiling_strategy_merge:
    """Merge all tiling strategies."""
    conda:
        "../envs/default.yaml"
    input:
        noTiles="results/tiling/{sample}/strategies/regions.noTiles.processed.bed.gz",
        twoTiles="results/tiling/{sample}/strategies/regions.twoTiles.processed.bed.gz",
        centeredTiles="results/tiling/{sample}/strategies/regions.centeredTiles.processed.bed.gz",
    output:
        "results/tiling/{sample}/regions.tiles.bed.gz",
    log:
        "results/logs/tiling/strategy_merge.{sample}.log",
    shell:
        """
        cat {input.noTiles} {input.twoTiles} {input.centeredTiles}  | \
        bgzip -dc  | \
        sort -k1,1 -k2,2n | \
        bgzip -c > {output} 2> {log}
        """
