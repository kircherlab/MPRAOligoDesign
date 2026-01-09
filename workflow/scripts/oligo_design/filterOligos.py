#!/usr/bin/env python

import pysam
import pandas as pd
import click
from filter import seqs_filter, regions_filter


@click.command()
@click.option(
    "--seqs",
    "seqs_file",
    required=False,
    type=click.Path(exists=True, readable=True),
    help="File containing the designed oligos ",
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
    default=False,
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
):
    repeatIndex = pysam.TabixFile(simple_repeats_file)
    TSSIndex = pysam.TabixFile(tss_pos_file)
    CTCFIndex = pysam.TabixFile(ctcf_motif_file)

    fail_reasons = {}
    failed = {}

    if regions_file:
        failed["regions"], fail_reasons = regions_filter(
            regions_file, repeatIndex, TSSIndex, CTCFIndex, repeat
        )

    if seqs_file:
        failed["seqs"], reasons = seqs_filter(seqs_file, maxHomLength)
        fail_reasons |= reasons

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
    seq_map = pd.read_csv(map_file, sep="\t")

    total = len(
        pd.concat([variant_map["REF_ID"], variant_map["ALT_ID"], seq_map["ID"]]).unique()
    )

    if "regions" in failed:
        seq_map = seq_map.drop(seq_map[seq_map["Region"].isin(failed["regions"])].index)
        variant_map = variant_map.drop(
            variant_map[variant_map["Region"].isin(failed["regions"])].index
        )

    if "seqs" in failed:
        seq_map = seq_map.drop(seq_map[seq_map["ID"].isin(failed["seqs"])].index).copy()

        variant_map = variant_map.drop(
            variant_map[variant_map["REF_ID"].isin(failed["seqs"])].index
        )
        variant_map = variant_map.drop(
            variant_map[variant_map["ALT_ID"].isin(failed["seqs"])].index
        )
        # remove alt ids if the associated ref id fails:
        ref_alt_dict = variant_map.groupby('REF_ID')['ALT_ID'].apply(list).to_dict()
        failed_ref_alt_seqs = flatten_list([ref_alt_dict[name] for name in failed["seqs"] if name in ref_alt_dict.keys()])
        seq_map = seq_map.drop(seq_map[seq_map["ID"].isin(failed_ref_alt_seqs)].index).copy()

        if remove_regions_without_variants:
            seq_map = seq_map.drop(
                seq_map[
                    ~seq_map["ID"].isin(
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

    seq_map.to_csv(map_out_file, compression="gzip", sep="\t", index=False)

    removed = total - len(
        pd.concat([variant_map["REF_ID"], variant_map["ALT_ID"], seq_map["ID"]]).unique()
    )

    return (total, removed)


if __name__ == "__main__":
    cli()
