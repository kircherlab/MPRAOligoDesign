#!/usr/bin/env python

import pysam
import pandas as pd
import click
import gzip
from filter import seqs_filter, regions_filter


@click.command()
@click.option(
    "--seqs",
    "seqs_file",
    required=False,
    type=click.Path(exists=True, readable=True),
    help="File containing the designed oligos",
)
@click.option(
    "--regions",
    "regions_file",
    required=False,
    type=click.Path(exists=True, readable=True),
    help="File containing the regions",
)
@click.option(
    "--variant-map",
    "variant_map_in",
    required=False,
    type=click.Path(exists=True, readable=True),
    help="Map that maps variants to region to ref and alt sequences",
)
@click.option(
    "--map",
    "map_in",
    required=False,
    type=click.Path(exists=True, readable=True),
    help="Map that maps sequences to regions",
)
@click.option(
    "--remove-regions-without-variants/--keep-regions-without-variants",
    "remove_regions_without_variants",
    required=False,
    help="Remove regions without variants",
)
@click.option(
    "--max_homopolymer_length",
    "maxHomLength",
    default=10,
    type=int,
    help="Maximum homopolymer length (def 10)",
)
@click.option(
    "--repeat",
    "repeat",
    default=0.25,
    type=float,
    help="Maximum fraction explained by a single simple repeat annotation",
)
@click.option(
    "--simple-repeats",
    "simple_repeats_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Bedfile for simple repeats",
)
@click.option(
    "--tss-positions",
    "tss_pos_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Bedfile for TSS positions",
)
@click.option(
    "--ctcf-motifs",
    "ctcf_motif_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Bedfile for ctcf motifs",
)
@click.option(
    "--output-map",
    "map_out",
    required=True,
    type=click.Path(writable=True),
    help="Output file of filtered sequence map",
)
@click.option(
    "--output-variant-map",
    "variant_map_out",
    required=False,
    type=click.Path(writable=True),
    help="Output file of filtered variant to region to sequence map",
)
# TODO: accept output files for each of the filters (optional but check if the paths are valid)
@click.option(
    "--output-failed-ctcf-bed",
    "failed_ctcf_bed",
    required=False,
    type=click.Path(writable=True),
    help="Output file of failed CTCF sites",
)
@click.option(
    "--output-failed-homopolymer-bed",
    "failed_homopolymer_bed",
    required=False,
    type=click.Path(writable=True),
    help="Output file of failed homopolymer sequences",
)
@click.option(
    "--output-failed-simpleRepeats-bed",
    "failed_simpleRepeats_bed",
    required=False,
    type=click.Path(writable=True),
    help="Output file of failed simple repeats",
)
@click.option(
    "--output-failed-TSS-bed",
    "failed_TSS_bed",
    required=False,
    type=click.Path(writable=True),
    help="Output file of failed TSS sites",
)
@click.option(
    "--output-failed-restriction-site-bed",
    "failed_restriction_site_bed",
    required=False,
    type=click.Path(writable=True),
    help="Output file of sequences with a restriction site",
)
def cli(
    seqs_file,
    regions_file,
    variant_map_in,
    map_in,
    remove_regions_without_variants,
    maxHomLength,
    repeat,
    simple_repeats_file,
    tss_pos_file,
    ctcf_motif_file,
    variant_map_out,
    map_out,
    failed_ctcf_bed,
    failed_homopolymer_bed,
    failed_simpleRepeats_bed,
    failed_TSS_bed,
    failed_restriction_site_bed,
):
    repeatIndex = pysam.Tabixfile(simple_repeats_file)
    TSSIndex = pysam.Tabixfile(tss_pos_file)
    CTCFIndex = pysam.Tabixfile(ctcf_motif_file)

    fail_reasons = {}
    failed = {}

    if regions_file:
        failed_tuples, fail_reasons = regions_filter(
            regions_file, repeatIndex, TSSIndex, CTCFIndex, repeat
        )
        failed["regions"] = [x[0] for x in failed_tuples]
        # write failed regions to bed files
        if failed_ctcf_bed:
            fail_reason = "CTCF"
            write_failed_regions(failed_tuples, failed_ctcf_bed, regions_file, fail_reason)

        if failed_simpleRepeats_bed:
            fail_reason = "repeats"
            write_failed_regions(failed_tuples, failed_simpleRepeats_bed, regions_file, fail_reason)

        if failed_TSS_bed:
            fail_reason = "TSS"
            write_failed_regions(failed_tuples, failed_TSS_bed, regions_file, fail_reason)


    if seqs_file:
        failed_tuples, reasons = seqs_filter(seqs_file, maxHomLength)
        fail_reasons |= reasons
        failed["seqs"] = [x[0] for x in failed_tuples]
        # write failed sequences to bed files
        if failed_homopolymer_bed:
            reason = "homopolymer"
            write_failed_sequence_regions(failed_tuples,
                                          region_sequence_map_file_path=map_in,
                                          bed_file_path=failed_homopolymer_bed,
                                          regions_file=regions_file, fail_reason=reason,
                                          ID_col_name="ID", region_col_name="Region")

        # write failed sequences to bed files
        if failed_restriction_site_bed:
            reason = "restriction"
            write_failed_sequence_regions(failed_tuples,
                                          region_sequence_map_file_path=map_in,
                                          bed_file_path=failed_restriction_site_bed,
                                          regions_file=regions_file, fail_reason=reason,
                                          ID_col_name="ID", region_col_name="Region")

    total, removed = write_output(
        failed,
        variant_map_in,
        map_in,
        remove_regions_without_variants,
        variant_map_out,
        map_out,
    )

    print(
        """Failed sequences:
    %d due to homopolymers
    %d due to simple repeats
    %d due to TSS site overlap
    %d due to CTCF site overlap
    %d due to EcoRI or SbfI restriction site overlap"""
        % (
            fail_reasons.get("hompol", 0),
            fail_reasons.get("repeats", 0),
            fail_reasons.get("TSS", 0),
            fail_reasons.get("CTCF", 0),
            fail_reasons.get("restrictions", 0),
        )
    )
    print("Total failed: %d" % (removed))
    print("Total passed sequences: %d" % (total - removed))


def write_failed_regions(failed_tuples, bed_file_path, regions_file, fail_reason):
    """Write the failed region bed file for the given fail reason"""
    # get the failed regions
    failed_region_ids = [x[0] for x in failed_tuples if x[1] == fail_reason]
    if len(failed_region_ids) == 0: # if no regions failed, nothing needs to be written
        return True

    if bed_file_path is None:
        raise ValueError("No output file path given")

    with open(bed_file_path, "w") as f:
        regions = gzip.open(regions_file, 'rt')
        for region in regions:
            region_split = region.strip().split("\t")
            rchrom, rstart, rend, rid = region_split[0:4]
            if rid in failed_region_ids:
                f.write(region)

    if bed_file_path.endswith(".gz"):
        written_regions = pd.read_csv(bed_file_path, sep="\t", compression=None)
        written_regions.to_csv(bed_file_path, sep="\t", index=False)
    return True


def write_failed_sequence_regions(failed_tuples, region_sequence_map_file_path, bed_file_path, regions_file, fail_reason, ID_col_name="ID", region_col_name="Region"):
    """Write the failed regions of sequences with homopolymers or restriction sites to a bed file"""
    map_df = pd.read_csv(region_sequence_map_file_path, sep="\t")
    failed_id = [x[0] for x in failed_tuples if x[1] == fail_reason]
    failed_map = map_df.loc[map_df[ID_col_name].isin(failed_id)].copy()
    failed_region_ids = failed_map[region_col_name].unique()
    if len(failed_region_ids) == 0:
        return True

    if bed_file_path is None:
        raise ValueError("No output file path given")


    with open(bed_file_path, "w") as f:
        # get all regions that failed and write them from the regions file
        regions = gzip.open(regions_file, 'rt')
        for region in regions:
            region_split = region.strip().split("\t")
            rchrom, rstart, rend, rid = region_split[0:4]
            if rid in failed_region_ids:
                f.write(region)

    if bed_file_path.endswith(".gz"):
        written_regions = pd.read_csv(bed_file_path, sep="\t", compression=None)
        written_regions.to_csv(bed_file_path, sep="\t", index=False)
    return True


def flatten_list(list_of_lists):
    return [item for sublist in list_of_lists for item in sublist]

def write_output(
    failed,
    variant_map_file,
    map_file,
    remove_regions_without_variants,
    variant_map_out_file,
    map_out_file,
):
    if variant_map_file:
        variant_map = pd.read_csv(variant_map_file, sep="\t")
    else:
        variant_map = pd.DataFrame(columns=["Variant", "Region", "REF_ID", "ALT_ID"])
    map = pd.read_csv(map_file, sep="\t")

    total = len(
        pd.concat([variant_map["REF_ID"], variant_map["ALT_ID"], map["ID"]]).unique()
    )

    if "regions" in failed:
        map = map.drop(map[map["Region"].isin(failed["regions"])].index)
        variant_map = variant_map.drop(
            variant_map[variant_map["Region"].isin(failed["regions"])].index
        )

    if "seqs" in failed:
        map = map.drop(map[map["ID"].isin(failed["seqs"])].index).copy()

        variant_map = variant_map.drop(
            variant_map[variant_map["REF_ID"].isin(failed["seqs"])].index
        )
        variant_map = variant_map.drop(
            variant_map[variant_map["ALT_ID"].isin(failed["seqs"])].index
        )
        # remove alt ids if the associated ref id fails:
        ref_alt_dict = variant_map.groupby('REF_ID')['ALT_ID'].apply(list).to_dict()
        failed_ref_alt_seqs = flatten_list([ref_alt_dict[name] for name in failed["seqs"] if name in ref_alt_dict.keys()])
        map = map.drop(map[map["ID"].isin(failed_ref_alt_seqs)].index).copy()

        if remove_regions_without_variants:
            map = map.drop(
                map[
                    ~map["ID"].isin(
                        pd.concat(
                            [variant_map["REF_ID"], variant_map["ALT_ID"]]
                        ).unique()
                    )
                ].index
            )

    if variant_map_out_file:
        variant_map.to_csv(
            variant_map_out_file, compression="gzip", sep="\t", index=False
        )

    map.to_csv(map_out_file, compression="gzip", sep="\t", index=False)

    removed = total - len(
        pd.concat([variant_map["REF_ID"], variant_map["ALT_ID"], map["ID"]]).unique()
    )

    return (total, removed)


if __name__ == "__main__":
    cli()
