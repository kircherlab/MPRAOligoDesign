import click
import gzip
import copy as cp
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
              help='Min overlap of two oligos')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output bed file')
def cli(input_file, oligo_length, min_overlap, output_file):

    if is_gz_file(input_file):
        file = gzip.open(input_file, 'rt')
    else:
        file = open(input_file, 'r')

    output = open(output_file, "w")
    for line in file:
        
        line_split = line.strip().split('\t')
        start, end = line_split[1:3]
        start, end = int(start), int(end)
        tiles = [tile for tile in tileCenterRegion(start, end, oligo_length, min_overlap)]
        n_tiles=len(tiles)
        n_tiles_half=int(n_tiles/2)

        for i, tile in enumerate(tiles):
            print_line = cp.deepcopy(line_split)
            print_line[1] = str(tile[0])
            print_line[2] = str(tile[1])

            # check for correct order. center must be in the center

            if i == 0:
                print_line[3] =  '%s_%s_tile%d-%d' % (print_line[3], STRAND_NAME[print_line[5]], n_tiles_half+1, n_tiles)
            elif (i> n_tiles_half):
                print_line[3] =  '%s_%s_tile%d-%d' % (print_line[3], STRAND_NAME[print_line[5]], i+1, n_tiles)
            else:
                print_line[3] =  '%s_%s_tile%d-%d' % (print_line[3], STRAND_NAME[print_line[5]], i, n_tiles)
            
            output.write('\t'.join(print_line) + '\n')
    output.close()

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

STRAND_NAME={'+': 'fwd', '-': 'rev', '.': 'none'}

def tileCenterRegion(start, end, oligo_length, min_overlap):
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
    for tile_left, tile_right in tileRegion(start, center_right, oligo_length, min_overlap, isLeft=True):
        yield (tile_left, tile_right)
    for tile_left, tile_right in tileRegion(center_left, end, oligo_length, min_overlap, isLeft=False):
        yield (tile_left, tile_right)


def tileRegion(start, end, oligo_length, min_overlap, isLeft=True):

    div = ((end - start) - min_overlap) // (oligo_length-min_overlap)
    if ((end - start) - min_overlap) % (oligo_length-min_overlap) != 0:
        div += 1

    # just center, left, right
    if div == 2:
        if isLeft:
            yield(start, start+oligo_length)
        else:
            yield(end-oligo_length, end)

    else:
        extension = (((oligo_length-min_overlap) - ((end - start) - min_overlap) %
                     (oligo_length-min_overlap))) % (oligo_length-min_overlap) // (div-1)
        first = ((oligo_length-min_overlap) - ((end - start) - min_overlap) % (oligo_length-min_overlap)) % (div-1)
        for i in range(div):
            tile_left = start
            tile_right = start + oligo_length
            # don't get center
            if not ((i == div-1 and isLeft) or (i == 0 and not isLeft)):
                yield(tile_left, tile_right)
            if (i == div-2 and isLeft) or (i == 0 and not isLeft):
                start = tile_right-(extension+min_overlap+first)
            else:
                start = tile_right-(extension+min_overlap)


if __name__ == '__main__':
    cli()
