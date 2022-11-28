def tiling_getInputFile(sample):
    """Return the input file for the given sample where tiling should be used."""
    if config["tiling"]["remove_edge_variants"]:
        return samples[sample]["bed_file"]
    else:
        return "results/tiling/%s/regions/extended.bed.gz" % sample
