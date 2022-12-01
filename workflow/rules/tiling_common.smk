def tiling_getInputFile(sample):
    """Return the input file for the given sample where tiling should be used."""
    if config["tiling"]["remove_edge_variants"]:
        return samples.loc[sample]["bed_file"]
    else:
        return "results/tiling/%s/regions/extended.bed.gz" % sample


def tiling_getMinTwoTiles():
    variant_edge_size = (
        config["tiling"]["variant_edge_exclusion"] * 2
        if config["tiling"]["strategies"]["two_tiles"]["include_variant_edge"]
        else 0
    )
    return min(
        [
            config["tiling"]["strategies"]["two_tiles"]["max"],
            2 * config["oligo_length"]
            - 2 * config["tiling"]["min_overlap"]
            - variant_edge_size,
        ]
    )
