#!/usr/bin/env python

import argparse

from types import *

import os
import re

from BCBio import GFF


def return_map(prefix, gene, transcript):

    if re.search('\-[0-9]+$', transcript):

        t_id = transcript.split('-')

        del t_id[-1]

    else:

        t_id = transcript

    a_id = prefix + '|' + ''.join(t_id)

    print('\t'.join(map(str, [gene, ''.join(t_id), a_id])))


def main():
    desc = """Generate a reference-based gene-transcript map across multiple CAT GFFs."""

    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(metavar='reference', type=str, nargs=1, dest='ref', help='CAT input GFF.')

    parser.add_argument(metavar='query', type=str, nargs='*', dest='qry', help='CAT output GFF.')

    parser.add_argument('--biotype', type=str, dest='biotype', action='append', default=['protein_coding', 'unknown_likely_coding'], help='Biotype to limit the map to.')

    parser.add_argument('--prefix', type=str, dest='prefix', action='append', help='Biotype to limit the map to.') # CHange

    args = parser.parse_args()

    print('### Starting ' + args.ref[0])

    ref_fh = open(args.ref[0], 'r')

    ref_name = os.path.splitext(args.ref[0])[0]

    limit_info = dict(gff_type=['gene', 'mRNA', 'CDS'])

    ref_gene_map = dict()

    for rec in GFF.parse(ref_fh, limit_info=limit_info):

        for gene in rec.features:

            ref_gene_biotype = gene.qualifiers.get('gene_biotype')

            if ref_gene_biotype[0] in args.biotype:

                ref_gene_locus = gene.qualifiers.get('locus_tag')

                ref_gene_map[gene.id] = ref_gene_locus[0]

                for transcript in gene.sub_features:

                    ref_transcript_locus = transcript.qualifiers.get('Name')

                    return_map(ref_name, ref_gene_locus[0], ref_transcript_locus[0])

                    del ref_transcript_locus

                del ref_gene_locus

    ref_fh.close()

    limit_info = dict(gff_type=['gene', 'transcript'])

    for qry in args.qry:

        qry_name = os.path.splitext(qry)[0]

        qry_fh = open(qry, 'r')

        print('### Starting ' + qry)

        for rec in GFF.parse(qry_fh, limit_info=limit_info):

            for gene in rec.features:

                qry_gene_biotype = gene.qualifiers.get('gene_biotype')

                if qry_gene_biotype[0] in args.biotype:

                    for transcript in gene.sub_features:

                        if transcript.qualifiers.get('transcript_class')[0] == 'ortholog' or transcript.qualifiers.get('transcript_class')[0] == 'putative_novel_isoform':

                            # Case 1 - Ortholog
                            # Solution: Get gene from current gene name and reference gene map
                            # Qualifier: transcript_class = ortholog
                            # Qualifier: transcript_class = putative_novel_isoform

                            qry_gene_id = transcript.qualifiers.get('source_gene')

                            qry_gene_locus = ref_gene_map[qry_gene_id[0]]

                            return_map(qry_name, qry_gene_locus, transcript.id)

                        elif transcript.qualifiers.get('transcript_class')[0] in ['poor_alignment', 'possible_paralog', 'putative_novel']:

                            # Case 2 - !Ortholog
                            # Solution: Gene gene from gene parent
                            # Qualifier: # transcript_class = poor_alignment
                            # Qualifier: # transcript_class = possible_paralog
                            # Qualifier: # transcript_class = putative_novel

                            qry_gene_locus = transcript.qualifiers.get('Parent')

                            return_map(qry_name, qry_gene_locus[0], transcript.id)

                        else:

                            print(transcript.id)


        qry_fh.close()


main()
