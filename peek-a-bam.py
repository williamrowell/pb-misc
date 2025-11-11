#!/usr/bin/env python3
# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "pysam",
# ]
# ///
"""
Check a PacBio BAM file for the presence of specific tags to infer the history.
"""

__version__ = '0.5.0'

import json
import math
import statistics
import pysam
import argparse

from collections import Counter
from typing import Any, Sequence, cast


MAX_QV = 60  # maximum meaningful Phred-scaled quality value
KINETICS_TAGS = {
  'fi',
  'ri',
  'fp',
  'rp',  # consensus fwd and rev ipd and pw
  'ip',
  'pw',  # subread ipd and pw
}
BASE_MODIFICATION_TAGS = {
  'MM',
  'ML',
  'Mm',
  'Ml',  # pre-standardization
}
SEGMENT_READS_TAGS = {
  'di',
  'dl',
  'dr',
  'ds',
}  # index, left and right adapater, and binary blob
DEMULTIPLEXED_TAGS = {
  'bc',
  'bq',  # barcode calls and quality
  'bl',
  'bt',  # leading/trailing barcode sequence
  'ql',
  'qt',  # leading/trailing barcode quality
  'bx',
  'ls',  # sequence lengths, binary blob of clipped sequences
}
HAPLOTAG_TAGS = {'HP', 'PS', 'PC'}  # haplotype tag, phase set, and phase confidence


def get_min_kmer(sequence: str) -> str:
  """Get the lexicographically smallest k-mer from a sequence length k."""
  rev_comp = sequence[::-1].translate(str.maketrans('ATCG', 'TAGC'))
  return min(sequence, rev_comp)


def rq_from_bq(bq: Sequence[Any]) -> int:
  """Compute read quality from an array of base qualities; cap at Q60"""
  read_len = len(bq)
  expectedErrors = sum([math.pow(10, -0.1 * x) for x in bq])
  try:
    return min(MAX_QV, math.floor(-10 * math.log10(expectedErrors / read_len)))
  except ZeroDivisionError:
    return -1


def error_rate_from_rqv(rqv: float) -> float:
  """Compute error rate from Phred-scaled read quality"""
  return math.pow(10, -0.1 * rqv) if rqv >= 0 else 1.0


def rqv_from_error_rate(error_rate: float) -> float:
  """Compute Phred-scaled read quality from error rate"""
  return -10 * math.log10(error_rate) if error_rate < 1.0 else -1.0


def parse_bam_header(bam_file: pysam.AlignmentFile) -> dict:
  """Parse the BAM header and return a dictionary."""
  header = cast(dict[str, Any], bam_file.header)
  read_groups = sorted(header.get('RG', []), key=lambda rg: rg.get('ID', ''))
  return {
    'is_aligned': bool(bam_file.nreferences),
    'read_groups': read_groups,
    'read_groups_ids': list({rg.get('ID', 'N/A') for rg in read_groups}),
    'movies': {rg.get('ID', 'N/A'): rg.get('PU', 'N/A') for rg in read_groups},
    'instrument_model': {rg.get('ID', 'N/A'): rg.get('PM', 'N/A') for rg in read_groups},
    'libraries': {rg.get('ID', 'N/A'): rg.get('LB', 'N/A') for rg in read_groups},
    'samples': {rg.get('ID', 'N/A'): rg.get('SM', 'N/A') for rg in read_groups},
  }


def calculate_movie_length(max_frame: int, framerate: int) -> str:
  """Calculate movie length in hours from max frame and framerate."""
  if max_frame == 0:
    return '0h'
  else:
    return f'{round((max_frame / framerate) / (60 * 60))}h'


def check_bam_file(
  bam_file_path: str,
  n_records: int,
  barcode_length: int,
  max_kmer_freq: float,
  framerate: int,
) -> dict:
  """
  Check the first n_records of a BAM file for the presence of specific tags
  and infer the history of the BAM file.
  """
  # load the BAM file
  save = pysam.set_verbosity(0)  # suppress [E::idx_find_and_load]
  bam_file = pysam.AlignmentFile(bam_file_path, 'rb', check_sq=False)
  pysam.set_verbosity(save)  # restore warnings

  header_info = parse_bam_header(bam_file)

  # initialize a flag to check if the BAM file contains non-CCS reads
  ccs = True
  hifi = True
  by_strand = False

  # initialize sets to collect per-read tags and values
  max_frame = {rg.get('ID', 'N/A'): 0 for rg in header_info['read_groups']}
  unique_tags: set[str] = set()
  unique_bq_scores: set[int] = set()
  np_list: list[int] = list()
  ec_list: list[float] = list()
  mg_list: list[float] = list()
  readqv_list: list[float] = list()
  readlength_list: list[int] = list()

  # initialize counters for first and last bases, and records per read group
  first_last_bases: Counter = Counter()
  records_per_rg: Counter = Counter()

  # read at most the first n_records records from the BAM file
  for record_index, record in enumerate(bam_file):
    if record_index >= n_records:
      break

    records_per_rg[record.get_tag('RG')] += 1

    # ccs reads have names like <movie_name>/<zmw>/ccs[/(fwd|rev)]
    qname = record.query_name
    if not qname or 'ccs' not in qname:
      ccs = False
    # by-strand ccs adds 'fwd' or 'rev' to read names
    if qname and qname.endswith(('fwd', 'rev')):
      by_strand = True

    # add the phred-scaled read quality to the list of read qualities
    # use the "rq" tag if available, otherwise compute from base qualities
    if record.has_tag('rq'):
      errorrate = 1.0 - float(record.get_tag('rq'))
      readqv = MAX_QV if errorrate == 0 else math.floor(-10 * math.log10(errorrate))
    else:
      readqv = rq_from_bq(record.query_qualities if record.query_qualities is not None else [])
    readqv_list.append(readqv)

    # hifi reads are ccs reads with read quality >= 20
    if readqv < 20:
      hifi = False

    # the we tag is frame number for last base contributing to consensus
    # highest frame number for each read group used to estimate movie length in seconds
    if record.has_tag('we'):
      id = record.get_tag('RG')
      max_frame[id] = max(max_frame[id], int(record.get_tag('we')))

    # track the frequency of the first and last barcode_length bases from each read
    # this can be used to infer if the BAM file is potentially multiplexed
    seq = record.query_sequence
    if seq:
      first_last_bases[get_min_kmer(seq[:barcode_length])] += 1
      first_last_bases[get_min_kmer(seq[-barcode_length:])] += 1

    # track the read lengths, base qualities, and tags
    readlength_list.append(record.query_length)
    if record.query_qualities is not None:
      unique_bq_scores.update(list(record.query_qualities))

    # unique tags are compared to known tags to infer record characteristics
    tags = getattr(record, 'tags', [])
    unique_tags.update(tag[0] for tag in tags)

    # for ccs reads, track the number of passes and effective coverage
    if record.has_tag('np'):
      np_list.append(int(record.get_tag('np')))
    if record.has_tag('ec'):
      ec_list.append(float(record.get_tag('ec')))

    # for aligned reads, track the gap-compressed identity
    if record.has_tag('mg'):
      mg_list.append(float(record.get_tag('mg')))

  output: dict[str, Any] = dict()

  # report on the movies, instrument models, libraries, and samples
  output['num_records_checked'] = record_index
  output['movies'] = header_info['movies']
  output['instrument_model'] = header_info['instrument_model']
  output['libraries'] = header_info['libraries']
  output['samples'] = header_info['samples']
  output['movie_length'] = dict(
    sorted(
      {id: calculate_movie_length(max_frame[id], framerate) for id in header_info['read_groups_ids']}.items(),
      key=lambda x: x[0],
    )
  )
  output['records_per_rg'] = dict(
    sorted(
      {id: f'{count} reads, {(count / record_index):.2%}' for id, count in records_per_rg.items()}.items(),
      key=lambda x: x[0],
    )
  )

  # report on the data type
  output['is_ccs'] = ccs
  output['is_by_strand'] = by_strand
  output['is_hifi'] = hifi
  output['is_segment_reads'] = bool(unique_tags & SEGMENT_READS_TAGS)

  # check for evidence of mux/demux
  output['is_demultiplexed'] = bool(unique_tags & DEMULTIPLEXED_TAGS)
  output['barcodes'] = list({rg.get('BC', 'N/A') for rg in header_info['read_groups']})
  # Does any sequence occur in the first/last barcode_length bases of more than
  # max_kmer_frequency of the reads? If so, it is potentially multiplexed.
  # It's also possible that this is an amplicon library.
  output['is_potentially_multiplexed'] = first_last_bases.most_common(1)[0][1] / n_records > max_kmer_freq
  if output['is_potentially_multiplexed']:
    output['common_barcodes'] = [f'{seq} ({count})' for seq, count in first_last_bases.most_common(10)]

  # check for evidence of alignment and haplotagging
  output['is_aligned'] = header_info['is_aligned']
  output['is_haplotagged'] = bool(unique_tags & HAPLOTAG_TAGS)

  # check tags for evidence of kinetics and base modification
  output['has_kinetics'] = bool(unique_tags & KINETICS_TAGS)
  output['has_base_modification'] = bool(unique_tags & BASE_MODIFICATION_TAGS)

  # check the number of bq bins; revio BAMs bin bq scores
  output['base_quality_bins'] = len(unique_bq_scores)

  # report read statistics
  error_rates = [error_rate_from_rqv(rqv) for rqv in readqv_list]
  output['read_quality'] = {
    'mean': round(rqv_from_error_rate(statistics.mean(error_rates)), 2),
    'min': round(rqv_from_error_rate(max(error_rates)), 2),
    'median': round(rqv_from_error_rate(statistics.median(error_rates)), 2),
    'max': round(rqv_from_error_rate(min(error_rates)), 2),
  }
  output['read_length'] = {
    'mean': round(statistics.mean(readlength_list), 2),
    'pstdev': round(statistics.pstdev(readlength_list), 2),
    'min': round(min(readlength_list)),
    'median': round(statistics.median(readlength_list), 1),
    'max': round(max(readlength_list)),
  }
  output['num_passes'] = (
    {
      'mean': round(statistics.mean(np_list), 2),
      'pstdev': round(statistics.pstdev(np_list), 2),
      'min': round(min(np_list)),
      'median': round(statistics.median(np_list), 1),
      'max': round(max(np_list)),
    }
    if np_list
    else None
  )
  output['effective_coverage'] = (
    {
      'mean': round(statistics.mean(ec_list), 2),
      'pstdev': round(statistics.pstdev(ec_list), 2),
      'min': round(min(ec_list)),
      'median': round(statistics.median(ec_list), 2),
      'max': round(max(ec_list)),
    }
    if ec_list
    else None
  )
  output['gap_compressed_identity'] = (
    {
      'mean': round(statistics.mean(mg_list), 2),
      'min': round(min(mg_list)),
      'median': round(statistics.median(mg_list), 2),
      'max': round(max(mg_list)),
    }
    if mg_list
    else None
  )

  # Close the BAM file
  bam_file.close()

  return output


def main(args: argparse.Namespace) -> None:
  output = check_bam_file(
    args.bam_file_path,
    args.n_records,
    args.barcode_length,
    args.max_kmer_freq,
    args.framerate,
  )
  print(json.dumps(output, indent=2))


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('bam_file_path', help='Path to the BAM file.')
  parser.add_argument(
    '--n_records',
    type=int,
    default=10000,
    help='Number of records to read from the BAM file',
  )
  parser.add_argument(
    '--barcode_length',
    type=int,
    default=16,
    help='Barcode length',
  )
  parser.add_argument(
    '--max_kmer_freq',
    type=float,
    default=0.05,
    help='Maximum frequency threshold for barcode check',
  )
  parser.add_argument(
    '--framerate',
    type=int,
    default=100,
    help='Framerate for movie length estimation',
  )
  parser.add_argument('--version', action='version', version=__version__)
  args = parser.parse_args()
  main(args)
