#!/usr/bin/env python3
import pytest
import importlib

mbc = importlib.import_module('mod_base_count')


def test_discretize_probability():
  assert mbc.discretize_probability(0.0) == 0
  assert mbc.discretize_probability(0.5) == 128
  assert mbc.discretize_probability(0.99999) == 255
  with pytest.raises(Exception) as e_info:  # noqa: F841
    mbc.discretize_probability(1.0)
    mbc.discretize_probability(-0.1)


def test_base_count():
  class MockRead:
    def __init__(self, seq):
      self.query_sequence = seq

  read = MockRead('ATCGATCGAA')
  assert mbc.base_count(read, 'A') == 4
  assert mbc.base_count(read, 'T') == 2
  assert mbc.base_count(read, 'C') == 2
  assert mbc.base_count(read, 'G') == 2
  assert mbc.base_count(read, 'g') == 2
  assert mbc.base_count(read, 'N') == 0


def test_modified_base_count():
  class MockRead:
    def __init__(self, mod_bases):
      self.modified_bases = mod_bases

  read = MockRead(
    {
      ('A', 0, 'a'): [(10, 200), (20, 100), (30, 255)],
      ('C', 0, 'm'): [(15, 150)],
    }
  )
  assert mbc.modified_base_count(read, ('A', 0, 'a'), 128) == 2
  assert mbc.modified_base_count(read, ('A', 0, 'a'), 200) == 2
  assert mbc.modified_base_count(read, ('A', 0, 'a'), 255) == 1
  assert mbc.modified_base_count(read, ('C', 0, 'm'), 128) == 1
  assert mbc.modified_base_count(read, ('C', 0, 'm'), 200) == 0
  assert mbc.modified_base_count(read, ('T', 1, 'a'), 128) == 0
