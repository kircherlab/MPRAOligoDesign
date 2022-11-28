import click

# options


@click.command()
@click.option('--input',
              'input_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Input bed file')
@click.option('--oligo-length',
              'oligo_length',
              required=True,
              type=int,
              help='Length of oligo')
@click.option('--min-overlap',
              'min_overlap',
              required=True,
              type=int,
              help='Length of oligo')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output bed file')
def cli(input_file, oligo_length, min_overlap, output_file):

    for line in open(input_file):
        chrom, start, end = line.strip().split('\t')
        start, end = int(start), int(end)
        for tile in tileRegion(start, end, oligo_length, min_overlap):
            print(chrom, tile[0], tile[1], sep='\t', file=open(output_file, 'a'))



def tileRegion(start, end, oligo_length, min_overlap):
    """Tiling function for a single region.

    Parameters
    ----------
    start : int
        Start of the region.
    end : int
        End of the region.
    oligo_length : int
        Length of the oligo.

    Returns
    -------
    generator
        touples of tiles
    """
    center = start + (end - start)/2
    center_left = int(center - (oligo_length/2))
    center_right = int(center + (oligo_length/2))

    yield (center_left, center_right)

    # left
    left_end = center_left + min_overlap

    div = (left_end - start) // (oligo_length-min_overlap/2)
    if (left_end - start) % (oligo_length-min_overlap/2):
        div += 1
    

    for i in range(start, end, oligo_length):
        yield (i, i + oligo_length)




if __name__ == '__main__':
    cli()
