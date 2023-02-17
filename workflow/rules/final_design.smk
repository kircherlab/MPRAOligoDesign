rule oligo_design_copy_variants_map:
    """Copy the variant region map to the final design folder."""
    conda:
        "../envs/default.yaml"
    input:
        "results/oligo_design/{sample}/design_variants_filtered.variant_region_map.tsv.gz",
    output:
        "results/final_design/{sample}/variant_region_map.tsv.gz",
    params:
        sample=lambda wc: wc.sample,
    log:
        "logs/final_design/copy_variants_map.{sample}.log",
    shell:
        """
        zcat {input} | \
        awk -v OFS='\\t' '{{if (NR == 1) {{print $0}}else{{print "{params.sample}:"$1,"{params.sample}:"$2,"{params.sample}:"$3,"{params.sample}:"$4}}}}' | \
        gzip -c > {output} 2> {log}
        """


rule oligo_design_copy_variants:
    """Copy the variants to the final design folder."""
    conda:
        "../envs/default.yaml"
    input:
        "results/oligo_design/{sample}/design_variants_filtered.variants.vcf.gz",
    output:
        "results/final_design/{sample}/variants.vcf.gz",
    params:
        sample=lambda wc: wc.sample,
    log:
        "logs/final_design/copy_variants.{sample}.log",
    shell:
        """
        zcat {input} | \
        awk -v OFS='\\t' '{{if ($0 ~ /^#/) {{print $0}}else{{$3="{params.sample}:"$3; print $0}}}}' | \
        bgzip -c > {output} 2> {log}
        """


rule final_design_copy_regions:
    """Copy the regions to the final design folder."""
    conda:
        "../envs/default.yaml"
    input:
        lambda wc: getFinalRegionFile(wc.sample),
    output:
        "results/final_design/{sample}/regions.bed.gz",
    params:
        sample=lambda wc: wc.sample,
    log:
        "logs/final_design/copy_regions.{sample}.log",
    shell:
        """
        zcat {input} | \
        awk -v OFS='\\t' '{{$4="{params.sample}:"$4;print $0}}' | \
        bgzip -c > {output} 2> {log}
        """


rule final_design_add_adapters:
    """
    Add adapters to the final designed sequences.
    """
    conda:
        "../envs/default.yaml"
    input:
        lambda wc: getFinalDesignFile(wc.sample),
    output:
        "results/final_design/{sample}/design.fa.gz",
    log:
        "logs/final_design/add_adapters.{sample}.log",
    params:
        left=config["oligo_design"]["adapters"]["left"],
        right=config["oligo_design"]["adapters"]["right"],
        sample=lambda wc: wc.sample,
    shell:
        """
        awk 'BEGIN{{
            seq="";header=""
        }}{{
            if ($0 ~ /^>/) {{
                if (seq != "") {{
                    print header;
                    print "{params.left}"seq"{params.right}";
                }}
                header="{params.sample}:"substr($0,2);
                seq="";
            }} else {{
                seq+=toupper($1);
            }}
        }}END{{
            print header;
            print "{params.left}"seq"{params.right}";
        }}' {input} | \ 
        bgzip -c > {output} 2> {log}
        """


rule final_design_combine_designs:
    """
    Combine designs
    """
    conda:
        "../envs/default.yaml"
    input:
        designs=expand(
            "results/final_design/{sample}/design.fa.gz", sample=datasets.index
        ),
    output:
        "results/final_design/design.fa.gz",
    log:
        "logs/final_design/combine_designs.log",
    shell:
        """
        zcat {input.designs} | bgzip -c > {output} 2> {log};
        """


rule final_design_combine_regions:
    """
    Combine designed regions
    """
    conda:
        "../envs/default.yaml"
    input:
        regions=expand(
            "results/final_design/{sample}/regions.bed.gz",
            sample=getRegionDatasets(),
        ),
    output:
        "results/final_design/regions.bed.gz",
    log:
        "logs/final_design/combine_regions.log",
    shell:
        """
        zcat {input.regions} | sort -k 1,1 -k2,2n | \
        bgzip -c > {output} 2> {log};
        """


rule final_design_combine_variants_idx:
    """
    Combine designed variants
    """
    conda:
        "../envs/default.yaml"
    input:
        "results/final_design/{sample}/variants.vcf.gz",
    output:
        "results/final_design/{sample}/variants.vcf.gz.tbi",
    log:
        "logs/final_design/combine_variants_idx.{sample}.log",
    shell:
        """
        tabix {input} &> {log};
        """


rule final_design_combine_variants:
    """
    Combine designed variants
    """
    conda:
        "../envs/default.yaml"
    input:
        variants=expand(
            "results/final_design/{sample}/variants.vcf.gz",
            sample=getVariantDatasets(),
        ),
        idx=expand(
            "results/final_design/{sample}/variants.vcf.gz.tbi",
            sample=getVariantDatasets(),
        ),
    output:
        "results/final_design/variants.vcf.gz",
    log:
        "logs/final_design/combine_variants.log",
    shell:
        """
        bcftools merge {input.variants} | \
        bgzip -c > {output} 2> {log};
        """


rule final_design_combine_variant_map:
    """
    Combine designed variant region map
    """
    conda:
        "../envs/default.yaml"
    input:
        variant_region_map=expand(
            "results/final_design/{sample}/variant_region_map.tsv.gz",
            sample=getVariantDatasets(),
        ),
    output:
        "results/final_design/variant_region_map.tsv.gz",
    log:
        "logs/final_design/combine_variant_map.log",
    shell:
        """
        (
            zcat {input.variant_region_map} | awk 'NR==1';
            for i in {input.variant_region_map}; do
                zcat $i | awk 'NR>1';
            done
        ) | gzip -c > {output} 2> {log};
        """
