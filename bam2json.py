#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "pysam",
# ]
# ///
"""
Parse BAM/SAM records into JSON structures.
"""

import argparse
import json
import sys
from typing import Optional
import pysam
import array as _array
from collections.abc import Sequence


__version__ = '0.1.0'


# SAM flag meanings
FLAG_MEANINGS = {
  0x1: 'read_paired',
  0x2: 'proper_pair',
  0x4: 'read_unmapped',
  0x8: 'mate_unmapped',
  0x10: 'read_reverse_strand',
  0x20: 'mate_reverse_strand',
  0x40: 'first_in_pair',
  0x80: 'second_in_pair',
  0x100: 'not_primary_alignment',
  0x200: 'read_fails_quality_checks',
  0x400: 'read_is_duplicate',
  0x800: 'supplementary_alignment',
}

# LocalContextFlags
LOCAL_CONTEXT_MEANINGS = {
  0x1: 'ADAPTER_BEFORE',
  0x2: 'ADAPTER_AFTER',
  0x4: 'BARCODE_BEFORE',
  0x8: 'BARCODE_AFTER',
  0x10: 'FORWARD_PASS',
  0x20: 'REVERSE_PASS',
  0x40: 'ADAPTER_BEFORE_BAD',
  0x80: 'ADAPTER_AFTER_BAD',
}


def infer_mode(path: str, fmt: Optional[str]) -> str:
  """
  Return pysam.AlignmentFile mode ('r' for SAM, 'rb' for BAM) based on
  explicit --input-format, filename extension, or fallback.
  """
  if fmt:
    if fmt.lower() == 'sam':
      return 'r'
    elif fmt.lower() == 'bam':
      return 'rb'
    else:
      raise ValueError("--input-format must be 'sam' or 'bam'")
  # If not provided, infer from extension when possible
  p = (path or '').lower()
  if p.endswith('.sam'):
    return 'r'
  if p.endswith('.bam'):
    return 'rb'
  # Default to BAM for unknown/STDIN; users can override with --input-format
  return 'rb'


def parse_flags(flag: int, flag_definitions: dict[int, str]) -> str:
  """
  Generic bitmask flag parser.

  Given an integer bitmask (`flag`) and a mapping of bit values to
  their string meanings (`flag_definitions`), return a comma-delimited
  string describing all the active flags.
  """
  meanings = [desc for bit, desc in flag_definitions.items() if flag & bit]
  return ','.join(meanings)


def _sanitize_tag_value(v):
  """
  Convert tag values into JSON-friendly Python types.

  If v is a (value, typecode) tuple as returned by pysam get_tag(..., with_value_type=True),
  use the SAM type code to coerce to an appropriate Python type:

    - 'A' or 'Z' -> str
    - integer types ('cCsSiI') -> int
    - 'f' -> float
    - 'B' -> array: return comma-delimited string of list of sanitized elements
    - 'H' -> hex byte-array: return hex string
    - None -> None

  If v is raw bytes/bytearray/memoryview, try decode utf-8, else return hex string.
  Lists/tuples are recursively sanitized.
  """
  # support pysam (value, typecode) tuple
  if isinstance(v, tuple) and len(v) == 2 and isinstance(v[1], str):
    val, tcode = v
    if val is None:
      return None
    t = tcode
    # string/char
    if t in ('A', 'Z'):
      return str(val)
    # integer-like types
    if t in ('c', 'C', 's', 'S', 'i', 'I'):
      return int(val)
    # float
    if t == 'f':
      return float(val)
    # general array: pysam returns (list, 'B') where list elements are primitives
    if t == 'B':
      # return a JSON-serializable list
      if val is None:
        return []
      # convert array.array and numpy-like sequences to plain lists
      if isinstance(val, _array.array):
        val = list(val)
      if hasattr(val, 'tolist') and not isinstance(val, (str, bytes, bytearray, memoryview)):
        try:
          val = val.tolist()
        except Exception:
          pass
      return [_sanitize_tag_value(x) for x in val] if val else []
    # hex byte array -> represent as hex string
    if t == 'H':
      if isinstance(val, (memoryview, bytes, bytearray)):
        return bytes(val).hex()
      return str(val)
    # fallback: sanitize underlying value
    return _sanitize_tag_value(val)

  # original handling for raw values
  if v is None:
    return None
  if isinstance(v, memoryview):
    v = v.tobytes()
  if isinstance(v, (bytes, bytearray)):
    try:
      s = bytes(v).decode('utf-8')
    except Exception:
      return bytes(v).hex()
    if s.lstrip('+-').isdigit():
      try:
        return int(s)
      except Exception:
        pass
    return s
  # convert array.array -> list, and numpy-like objects via tolist()
  if isinstance(v, _array.array):
    return [_sanitize_tag_value(x) for x in list(v)]
  if hasattr(v, 'tolist') and not isinstance(v, (str, bytes, bytearray, memoryview)):
    try:
      converted = v.tolist()
      return _sanitize_tag_value(converted)
    except Exception:
      pass
  if isinstance(v, (list, tuple, Sequence)):
    # Sequence includes lists/tuples; strings/bytes already handled above
    return [_sanitize_tag_value(x) for x in v]
  return v


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument(
    'input',
    nargs='?',
    default='-',
    help="Input SAM/BAM file path or '-' for stdin (default: '-')",
  )
  parser.add_argument(
    '-t',
    '--tags',
    metavar='TAG1,TAG2,...',
    help='Comma-delimited list of SAM optional tags to include; if omitted, include all tags present on each record.',
  )
  parser.add_argument(
    '--input-format',
    choices=['sam', 'bam'],
    help='Force input format (useful when reading from stdin). If omitted, '
    'the script attempts to infer from file extension, else defaults to BAM.',
  )
  args = parser.parse_args()

  tag_list = None
  if args.tags:
    tag_list = [t.strip() for t in args.tags.split(',') if t.strip()]

  mode = infer_mode(args.input, args.input_format)

  save = pysam.set_verbosity(0)
  # Open input (filename or stdin)
  try:
    af = pysam.AlignmentFile(args.input, mode, check_sq=False)
  except ValueError as e:
    # If we guessed BAM for stdin but it's actually SAM, try SAM as a fallback.
    if args.input == '-' and mode == 'rb':
      try:
        af = pysam.AlignmentFile(args.input, 'r')
      except Exception:
        raise e
    else:
      raise e
  pysam.set_verbosity(save)

  # Stream JSON array without holding everything in memory
  out = sys.stdout
  out.write('[')
  first = True

  try:
    for aln in af.fetch(until_eof=True):
      rec = {
        'query_name': aln.query_name,
        'flag': int(aln.flag),
        'flag_parsed': parse_flags(int(aln.flag), FLAG_MEANINGS),
        'query_length': int(aln.query_length or 0),
      }

      if tag_list:
        # Include only requested tags, if present on this record
        for t in tag_list:
          if aln.has_tag(t):
            # prefer typed form if available
            try:
              val = aln.get_tag(t, with_value_type=True)
            except TypeError:
              val = aln.get_tag(t)
            rec[t] = _sanitize_tag_value(val)

      else:
        # Include all tags. iterate tag names from aln.tags but fetch typed value
        for t, _ in aln.tags or []:
          try:
            val = aln.get_tag(t, with_value_type=True)
          except TypeError:
            val = aln.get_tag(t)
          rec[t] = _sanitize_tag_value(val)
      # Parse local context flag 'cx' if possible
      if 'cx' in rec:
        try:
          cx_val = rec['cx']
          if cx_val is not None:
            cx_int = int(cx_val)
            rec['cx_parsed'] = parse_flags(cx_int, LOCAL_CONTEXT_MEANINGS)
        except Exception:
          # ignore non-integer cx values
          pass

      if first:
        first = False
      else:
        out.write(',')
      json.dump(rec, out)
  finally:
    af.close()
    out.write(']\n')
    out.flush()


if __name__ == '__main__':
  main()
