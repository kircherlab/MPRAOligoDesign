import click
from Bio.Seq import Seq
import vcfpy
import pandas as pd
import pyranges as pr
import copy as cp
from pyfaidx import Fasta
# options


@click.command()
@click.option('--input-regions',
              'input_region_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Input region bed file')
@click.option('--input-variants',
              'input_variant_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Input variant vcf file')
@click.option('--reference',
              'reference_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Reference fasta file')
@click.option('--output-variants',
              'output_variant_file',
              required=True,
              type=click.Path(writable=True),
              help='Output variant file')
@click.option('--output-variants-removed',
              'output_removed_variant_file',
              required=True,
              type=click.Path(writable=True),
              help='Output removed variant file')
@click.option('--output-regions',
              'output_region_file',
              required=True,
              type=click.Path(writable=True),
              help='Output region file')
@click.option('--output-regions-without-variants',
              'output_removed_region_file',
              required=True,
              type=click.Path(writable=True),
              help='Output regions without variants')
@click.option('--output-design',
              'output_design_file',
              required=True,
              type=click.Path(writable=True),
              help='Output design file')
@click.option('--output-design-map',
              'output_design_map_file',
              required=True,
              type=click.Path(writable=True),
              help='Output design map file')
@click.option('--variant-edge-exclusion',
              'variant_edge_exclusion',
              required=True,
              type=int,
              help="Don't insert variants on the edge of the region")
@click.option('--remove-regions-without-variants/--keep-regions-without-variants',
              'remove_regions_without_variants',
              required=True,
              help="Remove regions without variants")
@click.option('--use-most-centered-region-for-variant/--use-all-regions-for-variant',
              'use_most_centered_region_for_variant',
              required=True,
              help="Only use one region where the variant is most centered on")
def cli(input_region_file, input_variant_file, reference_file, output_variant_file, output_removed_variant_file, output_region_file, output_removed_region_file, output_design_file, output_design_map_file, variant_edge_exclusion, remove_regions_without_variants, use_most_centered_region_for_variant):

    vcf_reader = vcfpy.Reader.from_path(input_variant_file)
    vcf_removed_writer = vcfpy.Writer.from_path(output_removed_variant_file, vcf_reader.header)
    vcf_header = vcf_reader.header
    vcf_header.add_info_line({"ID": "Region", "Number": ".", "Type": "String", "Description": "Matched regions"})
    vcf_header.add_info_line({"ID": "ALT_ID", "Number": ".", "Type": "String", "Description": "Corresponding alt ID(s)"})
    vcf_header.add_info_line({"ID": "REF_ID", "Number": ".", "Type": "String", "Description": "Corresponding REF ID(s)"})
    vcf_writer = vcfpy.Writer.from_path(output_variant_file, vcf_header)
    df = pd.read_csv(input_region_file, header=None, sep="\t")
    df.rename({0: "Chromosome", 1: "Start", 2: "End", 3: "Name", 4: "Score", 5: "Strand"}, axis=1, inplace=True)
    regions = pr.PyRanges(df)
    reference = Fasta(reference_file)

    regions_output = pd.DataFrame()
    ref_sequences = pd.DataFrame(columns="ID Sequence".split())
    alt_sequences = pd.DataFrame(columns="ID Sequence".split())
    design_map = pd.DataFrame(columns="Variant Region REF_ID ALT_ID".split())
    # take variant edges into account
    regions_variant_edges = cp.deepcopy(regions)
    regions_variant_edges.Start += variant_edge_exclusion
    regions_variant_edges.End -= variant_edge_exclusion
    for variant_record in vcf_reader:

        region_slice = regions_variant_edges[variant_record.CHROM, (variant_record.POS-1):variant_record.POS]

        if len(region_slice) == 0:
            vcf_removed_writer.write_record(variant_record)
            continue

        # revert variant edges
        region_slice.Start -= variant_edge_exclusion
        region_slice.End += variant_edge_exclusion

        # FIXME does not work with + and - strand counts them as one. Maybe use group by?
        if use_most_centered_region_for_variant and len(region_slice) > 1:
            distance_to_center = abs(variant_record.POS - 1 - region_slice.Start - region_slice.lengths())
            region_slice = region_slice[distance_to_center == distance_to_center.min()]
            region_slice = region_slice.head(1)

        regions_output = pd.concat([regions_output, region_slice.as_df()])

        ref_ids = []
        alt_ids = []
        region_names = []
        for region_name, alt_seq, ref_seq in getSequences(reference, region_slice, variant_record):
            if (len(variant_record.ID) == 0):
                variant_id = variant_record.CHROM + "-" + str(variant_record.POS) + "-" + variant_record.REF + "-" + variant_record.ALT[0].value
            else:
                variant_id = variant_record.ID[0]
            alt_id = "ALT_"+region_name + "_" + variant_id
            ref_id = "REF_"+region_name
            ref_sequences = pd.concat([ref_sequences, pd.DataFrame({"ID": [ref_id], "Sequence": [ref_seq]})])
            alt_sequences = pd.concat([alt_sequences, pd.DataFrame(
                {"ID": [alt_id], "Sequence": [alt_seq]})])
            design_map = pd.concat([design_map, pd.DataFrame(
                {"Variant": [variant_id], "Region": [region_name], "REF_ID": [ref_id], "ALT_ID": [alt_id]})])

            region_names.append(region_name)
            ref_ids.append(ref_id)
            alt_ids.append(alt_id)

        variant_record.INFO["Region"] = region_names
        variant_record.INFO["REF_ID"] = ref_ids
        variant_record.INFO["ALT_ID"] = alt_ids
        vcf_writer.write_record(variant_record)
    design_map.drop_duplicates(inplace=True, ignore_index=True)
    ref_sequences.drop_duplicates(inplace=True, ignore_index=True)
    alt_sequences.drop_duplicates(inplace=True, ignore_index=True)
    vcf_writer.close()
    vcf_removed_writer.close()

    design_map.to_csv(output_design_map_file, sep="\t", index=False, header=True)

    if (remove_regions_without_variants):
        regions_output.drop_duplicates(ignore_index=True).to_csv(output_region_file, sep="\t", header=False, index=False)

    else:
        regions.as_df().to_csv(output_region_file, sep="\t", header=False, index=False)

    removed_regions = pd.merge(regions.as_df(), regions_output, how="outer", indicator=True)
    removed_regions = removed_regions[removed_regions._merge == "left_only"]
    removed_regions.drop(columns="_merge", inplace=True)

    removed_regions.to_csv(output_removed_region_file, sep="\t", header=False, index=False)

    with open(output_design_file, "w") as design_file_writer:
        for index, row in ref_sequences.iterrows():
            design_file_writer.write(">"+row["ID"]+"\n")
            design_file_writer.write(row["Sequence"]+"\n")
        for index, row in alt_sequences.iterrows():
            design_file_writer.write(">"+row["ID"]+"\n")
            design_file_writer.write(row["Sequence"]+"\n")


def getSequences(reference, regions, variant_record):
    for region_idx, region in regions.as_df().iterrows():

        complement = region["Strand"] == "-"

        sequence = reference[region["Chromosome"]][region["Start"]:region["End"]]
        variant_position = variant_record.POS - sequence.start

        ref_seq = sequence.seq
        for i, alt in enumerate(variant_record.ALT):
            if i > 0:
                raise Exception("Only one ALT allele is supported for variant %s" % variant_record)
            if alt.type != "SNV":
                raise Exception("Only SNVs are supported for variant %s" % variant_record)
            alt_nuc = alt.value

            alt_seq = ref_seq[:variant_position] + alt_nuc + ref_seq[(variant_position+1):]
            if complement:
                alt_seq = str(Seq(alt_seq).reverse_complement())
                ref_seq = str(Seq(ref_seq).reverse_complement())
        yield (region["Name"], alt_seq, ref_seq)


if __name__ == '__main__':
    cli()
