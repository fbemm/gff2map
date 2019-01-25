#!/usr/bin/env python

import argparse

from types import *

import os
import re

from BCBio import GFF

def invert_dict(d):

    return dict([(v, k) for k, v in d.items()])


def return_map(prefix, gene, transcript):

    if re.search('\-[0-9]+$', transcript):

        t_id = transcript.split('-')

        del t_id[-1]

    else:

        t_id = transcript

    if re.search('^hgmDDN_', gene):

        t_id = gene + '.t1'

    a_id = prefix + '|' + transcript

    print('\t'.join(map(str, [gene, ''.join(t_id), a_id])))


def read_hgm(hgm, ddn):

    hgm_fh = open(hgm, 'r').read().splitlines()

    id_map = {}

    map = {}

    for line in hgm_fh:

        if re.match(r'^#', line):

            header = re.split(r'\s+', line)

            id_map.update({header[1]: header[2]})

        else:

            members = filter(None, re.split(r'\s+', re.sub(r'[()]*', '', line)))

            ids = []

            acc = []

            for member in members:

                member_data = re.split(r',', member)

                acc.append(member_data[0])

                ids.append(id_map[member_data[0]] + "|" + member_data[1])

            if len(set(acc)) is len(ids) and invert_dict(id_map)[ddn] not in acc:

                new = 'hgmDDN_' + str( len(map) + 1 )

                for transcript in ids:

                    map.update({transcript: new})

    return map


def main():

    desc = """Generate a reference-based gene-transcript map across multiple CAT GFFs."""

    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(metavar='reference', type=str, nargs=1, dest='ref', help='CAT input GFF.')

    parser.add_argument(metavar='query', type=str, nargs='*', dest='qry', help='CAT output GFF.')

    parser.add_argument('--biotype', type=str, dest='biotype', action='append', default=['protein_coding', 'unknown_likely_coding'], help='Biotype to limit the map to.')

    parser.add_argument('--hgm', type=str, dest='hgm', nargs='*', help='Optional hgm file for de novo loci.')

    parser.add_argument('--ddn', type=str, dest='ddn', nargs='*', help='Optional reference name to de novo filtering.')

    args = parser.parse_args()

    if args.hgm[0]:

        hgm_map = read_hgm(args.hgm[0],args.ddn[0])

    ref_fh = open(args.ref[0], 'r')

    ref_name = os.path.splitext(os.path.basename(args.ref[0]))[0]

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

        qry_name = os.path.splitext(os.path.basename(qry))[0]

        qry_fh = open(qry, 'r')

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

                        elif transcript.qualifiers.get('transcript_class')[0] in ['putative_novel']:

                            # Case 2 - Maybe Ortholog but not in Reference
                            # Solution: Gene gene from gene parent
                            # Qualifier: # transcript_class = putative_novel

                            full_id = qry_name + "|" + transcript.id

                            if full_id in hgm_map:

                                qry_gene_locus = hgm_map[full_id]

                                return_map(qry_name, qry_gene_locus, transcript.id)

                            else:

                                qry_gene_locus = transcript.qualifiers.get('Parent')

                                return_map(qry_name, qry_gene_locus[0], transcript.id)

                        elif transcript.qualifiers.get('transcript_class')[0] in ['poor_alignment', 'possible_paralog']:

                            # Case 3 - !Ortholog
                            # Solution: Gene gene from gene parent
                            # Qualifier: # transcript_class = poor_alignment
                            # Qualifier: # transcript_class = possible_paralog

                            qry_gene_locus = transcript.qualifiers.get('Parent')

                            return_map(qry_name, qry_gene_locus[0], transcript.id)

                        else:

                            print(transcript.id)

        qry_fh.close()


main()
