#!/usr/bin/env python3
import pytest
import importlib
import sys
import json
from pathlib import Path

b2j = importlib.import_module('bam2json')


def test_infer_mode_and_parse_flags():
  assert b2j.infer_mode('file.sam', None) == 'r'
  assert b2j.infer_mode('file.bam', None) == 'rb'
  assert b2j.infer_mode('-', None) == 'rb'
  assert b2j.infer_mode('any', 'sam') == 'r'
  assert b2j.infer_mode('any', 'bam') == 'rb'
  with pytest.raises(ValueError):
    b2j.infer_mode('any', 'unknown')

  mapping = {0x1: 'one', 0x2: 'two', 0x4: 'four'}
  assert b2j.parse_flags(0x1, mapping) == 'one'
  assert b2j.parse_flags(0x3, mapping) == 'one,two'
  assert b2j.parse_flags(0x0, mapping) == ''


class FakeAlignment:
  """
  tags param: list of (tag, value, optional_typecode).
  We store an internal map and expose .tags as sequence of (tag, None) so the
  production code can iterate names; get_tag supports with_value_type=True.
  """

  def __init__(self, query_name, flag, query_length, tags):
    self.query_name = query_name
    self.flag = flag
    self.query_length = query_length
    self._tag_map = {}
    self.tags = []
    for item in tags:
      if len(item) == 3:
        tag, val, tcode = item
      else:
        tag, val = item
        tcode = None
      self._tag_map[tag] = (val, tcode)
      # keep second element for unpacking in main (ignored)
      self.tags.append((tag, None))

  def has_tag(self, t):
    return t in self._tag_map

  def get_tag(self, t, with_value_type=False):
    val, tcode = self._tag_map[t]
    if with_value_type:
      # If no explicit tcode, attempt to infer a sensible SAM-like code
      if tcode is not None:
        return (val, tcode)
      # infer
      if isinstance(val, (bytes, bytearray, memoryview)):
        return (val, 'H')
      if isinstance(val, list):
        return (val, 'B')
      if isinstance(val, float):
        return (val, 'f')
      if isinstance(val, int):
        return (val, 'i')
      if val is None:
        return (None, 'Z')
      return (str(val), 'Z')
    return val


def FakeAFFactory(aligns):
  def factory(path, mode, check_sq=False):
    class FA:
      def __init__(self, aligns):
        self._aligns = aligns
        self.closed = False

      def fetch(self, until_eof=True):
        for a in self._aligns:
          yield a

      def close(self):
        self.closed = True

    return FA(aligns)

  return factory


def _install_fake_pysam(monkeypatch, alignments):
  monkeypatch.setattr(b2j.pysam, 'AlignmentFile', FakeAFFactory(alignments))
  monkeypatch.setattr(b2j.pysam, 'set_verbosity', lambda v: 0)


def test_main_streams_many_tag_types_and_parses_cx(monkeypatch, capsys):
  a = FakeAlignment(
    query_name='read_mixed',
    flag=0x10,
    query_length=100,
    tags=[
      ('A', 'x', 'A'),  # char
      ('B', 'longstring', 'Z'),  # string
      ('C', 42, 'i'),  # int
      ('D', 3.14, 'f'),  # float
      ('E', b'\xde\xad', 'H'),  # hex bytes
      ('F', bytearray(b'\xca\xfe'), 'H'),
      ('G', [1, 2, 3], 'B'),  # array of ints
      ('H', None, 'Z'),  # explicit None
      ('cx', 3, 'i'),  # cx as integer
    ],
  )
  _install_fake_pysam(monkeypatch, [a])

  monkeypatch.setattr(sys, 'argv', ['bam2json.py', '-'])

  b2j.main()

  out = capsys.readouterr().out
  data = json.loads(out)
  assert len(data) == 1
  rec = data[0]
  assert rec['A'] == 'x'
  assert rec['B'] == 'longstring'
  assert rec['C'] == 42
  assert abs(rec['D'] - 3.14) < 1e-9
  assert rec['E'] == 'dead'
  assert rec['F'] == 'cafe'
  assert rec['G'] == [1, 2, 3]
  assert rec['H'] is None
  assert rec['cx'] == 3
  assert rec['cx_parsed'] == 'ADAPTER_BEFORE,ADAPTER_AFTER'


def test_main_tags_filter_preserves_empty_values(monkeypatch, capsys):
  a = FakeAlignment(
    query_name='rX',
    flag=0,
    query_length=50,
    tags=[('NM', '', 'Z')],
  )
  _install_fake_pysam(monkeypatch, [a])

  monkeypatch.setattr(sys, 'argv', ['bam2json.py', '-', '-t', 'NM,YY'])
  b2j.main()

  out = capsys.readouterr().out
  data = json.loads(out)
  assert len(data) == 1
  rec = data[0]
  assert 'NM' in rec
  assert rec['NM'] == ''
  assert 'YY' not in rec


def test_integration_with_sample_bam(capsys):
  repo_root = Path(__file__).resolve().parents[1]
  bams = list(repo_root.rglob('*.bam'))
  if not bams:
    pytest.skip('no .bam file found in repository for integration test')
  sample = str(bams[0])

  sys_argv = sys.argv
  try:
    sys.argv = ['bam2json.py', sample]
    b2j.main()
  finally:
    sys.argv = sys_argv

  out = capsys.readouterr().out
  data = json.loads(out)
  assert isinstance(data, list)
  for rec in data:
    assert 'query_name' in rec
    assert 'flag' in rec
    assert 'flag_parsed' in rec
    assert 'query_length' in rec
