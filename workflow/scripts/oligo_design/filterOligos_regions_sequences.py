#!/usr/bin/env python

import pysam
import pandas as pd
import click
from filter import seqs_filter, regions_filter

@click.command()
@click.option('--seqs',
              'seqs_file',
              required=False,
              type=click.Path(exists=True, readable=True),
              help='File containing the designed oligos ')
@click.option('--regions',
              'regions_file',
              required=False,
              type=click.Path(exists=True, readable=True),
              help='File containing the regions')
@click.option('--map',
              'map_in',
              required=False,
              type=click.Path(exists=True, readable=True),
              help='Map that maps sequences to regions')
@click.option('--max_homopolymer_length',
              'maxHomLength',
              default=10,
              type=int,
              help='Maximum homopolymer length (def 10)')
@click.option('--repeat',
              'repeat',
              default=0.25,
              type=float,
              help='Maximum fraction explained by a single simple repeat annotation')
@click.option('--simple-repeats',
              'simple_repeats_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Bedfile for simple repeats')
@click.option('--tss-positions',
              'tss_pos_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Bedfile for TSS positions')
@click.option('--ctcf-motifs',
              'ctcf_motif_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Bedfile for ctcf motifs')
@click.option('--output-map',
              'map_out',
              required=True,
              type=click.Path(writable=True),
              help='Output file of filtered sequence map')
def cli(seqs_file, regions_file, map_in, maxHomLength, repeat, simple_repeats_file, tss_pos_file, ctcf_motif_file, map_out):

    repeatIndex = pysam.Tabixfile(simple_repeats_file)
    TSSIndex = pysam.Tabixfile(tss_pos_file)
    CTCFIndex = pysam.Tabixfile(ctcf_motif_file)

    fail_reasons = {}
    failed = {}

    if regions_file:
        failed["regions"], fail_reasons = regions_filter(regions_file, repeatIndex, TSSIndex, CTCFIndex, repeat)

    if seqs_file:
        failed["seqs"], reasons = seqs_filter(seqs_file, maxHomLength)
        fail_reasons |= reasons
    
    total, removed = write_output(failed, map_in, map_out)

    print("""Failed sequences: 
    %d due to homopolymers 
    %d due to simple repeats 
    %d due to TSS site overlap 
    %d due to CTCF site overlap 
    %d due to EcoRI or SbfI restriction site overlap""" % (fail_reasons.get("hompol", 0), fail_reasons.get("repeats", 0), fail_reasons.get("TSS", 0), fail_reasons.get("CTCF", 0), fail_reasons.get("restrictions", 0)))
    print("Total failed: %d" % (removed))
    print("Total passed sequences: %d" % (total-removed))


def write_output(failed, map_file, out_file):

    map = pd.read_csv(map_file, sep="\t")

    total = len(map["ID"].unique())

    if "regions" in failed:
        for rid in failed["regions"]:
            if rid in map["Region"].values:
                map = map.drop(map[map["Region"] == rid].index)

    if "seqs" in failed:
        for sid in failed["seqs"]:
            if sid in map["ID"].values:
                map = map.drop(map[map["ID"] == sid].index)
    
    map.to_csv(out_file, compression='gzip', sep="\t", index=False)

    removed = total -  len(map["ID"].unique())

    return (total, removed)

if __name__ == '__main__':
    cli()

