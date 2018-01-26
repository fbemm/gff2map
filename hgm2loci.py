#!/usr/bin/env python

import argparse

from types import *

import os
import re

from BCBio import GFF

desc = """Generate a HomGenMapping-based gene-transcript map for CAT's de novo genes."""

parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(metavar='hgm', type=str, nargs=1, dest='hgm', help='HomGenMapping homologs file.')

parser.add_argument(metavar='gff', type=str, nargs='*', dest='gff', help='CAT output GFF.')

parser.add_argument('--biotype', type=str, dest='biotype', action='append',
                    default=['protein_coding', 'unknown_likely_coding'],
                    help='Biotype to limit the map to.')

args = parser.parse_args()


