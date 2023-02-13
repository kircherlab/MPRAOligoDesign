import math


rule centering_create_regions:
    conda:
        "../envs/default.yaml"
    input:
        variants=lambda wc: datasets.loc[wc.sample]["vcf_file"],
        genome_file=config["reference"]["genome"],
    output:
        regions="results/centering/{sample}/regions.centered.bed.gz",
    params:
        l=math.ceil(config["oligo_length"] / 2),
        r=math.floor(config["oligo_length"] / 2) - 1,
    log:
        "logs/centering/create_regions.{sample}.log",
    shell:
        """
        (
            if [[ $(file -b --mime-type "{input.variants}") == 'application/gzip' || $(file -b --mime-type "{input.variants}") == 'application/x-gzip' ]]; then
                zcat {input.variants}
            else
                cat {input.variants}
            fi
        ) | \
        grep -v '^#' | \
        awk -v "OFS=\\t" '{{print $1,$2-1,$2,$3,".","."}}' | \
        bedtools slop -l {params.l} -r {params.r} -i - -g {input.genome_file} | \
        sort -k1,1 -k2,2n | bgzip -c > {output} 2> {log}
        """
