#!/usr/bin/env python

import argparse

from types import *
from itertools import *

import os
import re
import sys

import pandas as pd


def get_gene(transcripts, gene_map):

    transcripts_set = re.split(r',', transcripts)

    gene_list = [gene_map[x] for x in transcripts_set]

    gene_set = set(gene_list)

    if len(gene_set) == 1:

        return list(gene_set)[0]

    else:

        return False


def main():

    desc = """Generate a transmap-guided expression matrix across multiple RSEM outputs."""

    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(metavar='transmap', type=str, nargs=1, dest='map', help='Input transmap.')

    parser.add_argument(metavar='matrix', type=str, nargs=1, dest='out', help='Output matrix.')

    parser.add_argument(metavar='rsem', type=str, nargs='*', dest='exp', help='Input RSEM.')

    args = parser.parse_args()

    map_lines = open(args.map[0], 'r').read().splitlines()

    gene_map = dict()

    for line in map_lines:

        gene,isoform,allele = re.split(r'\t+', line)

        gene_map.update({allele: gene})

    mat = pd.DataFrame(index=list(set(gene_map.values())))

    for exp in args.exp:

        exp_name = os.path.splitext(os.path.basename(exp))[0]

        exp_line = open(exp, 'r').read().splitlines()

        for line in islice(exp_line, 1, None):

            data = re.split(r'\t+', line)

            gene = get_gene(data[1], gene_map)

            exp_count = data[4]

            mat.at[gene, exp_name] = exp_count

    mat.to_csv(sys.stdout, sep='\t', na_rep=0)

main()