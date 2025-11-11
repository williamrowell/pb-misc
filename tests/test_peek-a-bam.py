#!/usr/bin/env python3
import pytest
import importlib

pab = importlib.import_module('peek-a-bam')


class MockHeader:
  def __init__(self, nreferences=0, header=None):
    self.nreferences = nreferences
    self.header = header


def test_get_min_kmer() -> None:
  assert pab.get_min_kmer('ACGT') == 'ACGT'
  assert pab.get_min_kmer('GGTT') == 'AACC'


def test_rq_from_bq() -> None:
  assert pab.rq_from_bq([30, 30, 30, 30]) == 30
  assert pab.rq_from_bq([29, 30, 30, 30]) == 29
  assert pab.rq_from_bq([1000, 1000, 1000, 1000]) == pab.MAX_QV
  assert pab.rq_from_bq([]) == -1


def test_error_rate_from_rqv() -> None:
  assert pab.error_rate_from_rqv(30.0) == pytest.approx(0.001)
  assert pab.error_rate_from_rqv(29.9) == pytest.approx(0.0010232929922807535)
  assert pab.error_rate_from_rqv(20.0) == pytest.approx(0.01)
  assert pab.error_rate_from_rqv(-1.0) == pytest.approx(1.0)


def test_rqv_from_error_rate() -> None:
  assert pab.rqv_from_error_rate(0.001) == pytest.approx(30.0)
  assert pab.rqv_from_error_rate(0.01) == pytest.approx(20.0)
  assert pab.rqv_from_error_rate(1.0) == pytest.approx(-1.0)
  assert pab.rqv_from_error_rate(0.0011) == pytest.approx(29.5861)


def test_calculate_movie_length() -> None:
  assert pab.calculate_movie_length(360000, 100) == '1h'
  assert pab.calculate_movie_length(0, 100) == '0h'
  assert pab.calculate_movie_length(720000, 200) == '1h'


@pytest.fixture
def mock_header_aligned() -> MockHeader:
  return MockHeader(
    nreferences=5,
    header={
      'RG': [
        {
          'ID': 'rg1',
          'PU': 'movie1',
          'PM': 'Sequel',
          'LB': 'lib1',
          'SM': 'sample1',
        },
        {
          'ID': 'rg2',
          'PU': 'movie2',
          'PM': 'Revio',
          'LB': 'lib2',
          'SM': 'sample2',
        },
      ]
    },
  )


def test_parse_mock_header_aligned(mock_header_aligned: MockHeader) -> None:
  header_info = pab.parse_bam_header(mock_header_aligned)
  assert header_info['is_aligned'] is True
  assert header_info['movies'] == {'rg1': 'movie1', 'rg2': 'movie2'}
  assert header_info['instrument_model'] == {'rg1': 'Sequel', 'rg2': 'Revio'}
  assert header_info['libraries'] == {'rg1': 'lib1', 'rg2': 'lib2'}
  assert header_info['samples'] == {'rg1': 'sample1', 'rg2': 'sample2'}


@pytest.fixture
def mock_header_unaligned() -> MockHeader:
  return MockHeader(
    nreferences=0,
    header={
      'RG': [
        {
          'ID': 'rg1',
          'PU': 'movie1',
          'PM': 'Sequel',
          'LB': 'lib1',
          'SM': 'sample1',
        },
      ]
    },
  )


def test_parse_mock_header_unaligned(mock_header_unaligned: MockHeader) -> None:
  header_info = pab.parse_bam_header(mock_header_unaligned)
  assert header_info['is_aligned'] is False
