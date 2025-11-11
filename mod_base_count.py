#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "pysam",
# ]
# ///
"""Assess base modifications in a BAM file.
Output: <read_name> <read_length> <total_base_contexts> <modified_base_contexts>
"""

import sys
import argparse
import pysam


__version__ = '0.1.0'


MODIFICATIONS = {
  '6mA': [('A', 0, 'a'), ('T', 1, 'a')],  # 6mA is scored asymmetrically
  '5mC': [('C', 0, 'm')],
}


def discretize_probability(prob: float) -> int:
  """Convert a probability in [0.0, 1.0) to a discrete integer in [0, 255]."""
  if not 0.0 <= prob < 1.0:
    raise ValueError('Probability must be in the range [0.0, 1.0).')
  return int(prob * 256)


def base_count(read: pysam.AlignedSegment, base: str) -> int:
  """Count the number of occurrences of a base in a read."""
  if read.query_sequence is None:
    return 0
  else:
    return read.query_sequence.count(base.upper())


def modified_base_count(read: pysam.AlignedSegment, modification: tuple[str, int, str], min_prob: int) -> int:
  """Count the number of modified bases in a read above a probability threshold."""
  if read.modified_bases and modification in read.modified_bases:
    return len([_ for _ in read.modified_bases[modification] if _[1] >= min_prob])
  else:
    return 0


def main():
  parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('bam', help='Path to BAM file')
  parser.add_argument(
    '-m',
    '--modification',
    choices=MODIFICATIONS.keys(),
    default='6mA',
    help='Type of base modification to assess.',
  )
  probability = parser.add_argument_group(
    'Score threshold (choose one):', 'Minimum score to consider a base as modified.'
  )
  exclusive = probability.add_mutually_exclusive_group(required=True)
  exclusive.add_argument('-p', '--min_prob', type=float, help='Minimum probability for base modification (e.g., 0.5).')
  exclusive.add_argument(
    '-i', '--min_int', type=int, help='Minimum integer representation for base modification (e.g., 127).'
  )
  parser.add_argument('--version', action='version', version=__version__)
  args = parser.parse_args()

  if args.min_prob:
    if args.min_prob and not 0.0 <= args.min_prob < 1.0:
      print('min_prob must be in the range [0.0, 1.0).', file=sys.stderr)
      sys.exit(1)
    min_prob = discretize_probability(args.min_prob)
  else:
    if not 0 <= args.min_int <= 255:
      print('min_int must be in the range [0, 255].', file=sys.stderr)
      sys.exit(1)
    min_prob = args.min_int
  bam = pysam.AlignmentFile(args.bam, check_sq=False)
  for read in bam:
    total_bases = 0
    mod_bases = 0
    for mod in MODIFICATIONS[args.modification]:
      base = mod[0]
      total_bases += base_count(read, base)
      mod_bases += modified_base_count(read, mod, min_prob)
    print('\t'.join([read.query_name, str(read.query_length), str(total_bases), str(mod_bases)]))
  bam.close()


if __name__ == '__main__':
  main()
