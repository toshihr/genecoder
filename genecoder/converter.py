# -*- coding: utf-8 -*-
'''Converter from IARC-TP53DATABASE to genecoder input file.
'''
from __future__ import absolute_import, division, print_function, unicode_literals
try:
    from future_builtins import zip
except ImportError:
    # Python 3 raise ImportError
    pass
import csv
import re
import sys
import itertools
import os

# old manually created format
L2L3 = dict(
    mode_use_mutation_column=False,
    filter_wildtype={'domain': ('Wild Type',)},
    filter_only_for={},
    filter_omit={},
    columns_output_left=('PubMed', 'mutation_id',),
    columns_output_right=('Codon_number', 'Sample_source', 'Tumor_origin', 'Morphology',
                          'Histological Grade', 'Stage', 'TNM', 'p53_IHC', 'HER2(c-erbB-2)',
                          'ER status', 'PR status', 'Tumor subtype',
                          'OS(months)', 'OS(event)', 'RFS(months)', 'RFS(event)'),
    columns_output_generated = ('region_name', 'seq_category', 'seq_aa', 'seq_na'),
    columns_error=('mutation_id',),
)

IARC = dict(
    mode_use_mutation_column=True,
    filter_wildtype={},
    filter_only_for={'Short_topo': ('COLON',), 'Tumor_origin': ('primary',)},
    filter_omit={'Type': ('del', 'CC tandem', 'complex', 'ins', 'NA', 'tandem'), },
    columns_output_left=('Mutation_ID', 'Short_topo'),
    columns_output_right=('Sample_Name', 'Sample_source', 'TNM', 'Grade', 'Stage', 'Age',
                          'Ethnicity', 'Population', 'Country', 'Ref_ID', 'PubMed'),
    columns_output_generated = ('region_name', 'seq_category', 'seq_aa', 'seq_na'),
    columns_error=('Mutation_ID',),
)

# column_mutation_NA: 変異情報(遺伝子上)のある列名 (変異情報例：c.123G>A)
column_mutation_NA = 'c_description'
# column_mutation_AA: 変異情報(タンパク上)のある列名 (変異情報例：p.R175H)
column_mutation_AA = 'ProtDescription'
# column_codon_number: コドン番号(1スタート), mutation_AAと同じ番号
column_codon_number = 'Codon_number'
# column_wildtype_codon: 変異前コドンを表す塩基３文字のある列名
column_wildtype_codon = 'WT_codon'
# column_mutant_codon: 変異後コドンを表す塩基３文字のある列名
column_mutant_codon = 'Mutant_codon'
# regions: region name & it's range in amino acid sequence (index notation
# is same as the paper, ex. (1,5) == seq[0,5] in python)
regions = {'Transactivation': (1, 50),
           'SH3': (63, 97),
           'DNA binding': (102, 292),
           'Tetramerization': (323, 356),
           'Regulatory domain': (363, 393),
           'L1S': (112, 141),
           'L2': (163, 195),
           'L2inner': (164, 194),
           'L3': (236, 251),
           'L3inner': (237, 250),
           'H': (278, 286),
           'I': (13, 23),
           'II': (117, 142),
           'III': (171, 181),
           'IV': (234, 258),
           'V': (270, 286)}


# ================================================ MAIN ============================================
__ROOT_DIR = os.path.join(os.path.dirname(__file__), 'rawdata')
# codon table
codon_table = dict(zip((''.join(a_list) for a_list in itertools.product('TCAG', repeat=3)),
                       'FFLLSSSSYYXXCCXWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
# regex for mutation format (the column 'c_description' for na_regex, 'ProtDescription' for aa_regex
# ) with index of pos, before, after
mutation_format_na = (re.compile('c\.([0-9]+)(.)>(.)'), (0, 1, 2))
mutation_format_aa = (re.compile('p\.(.)([0-9]+)(.)'), (1, 0, 2))


def NA_to_AA(seq):
    ''' convert nucleotide acid sequence to amino acid sequece. '''
    AA = []
    for i in range(len(seq) // 3):
        codon = seq[i * 3:(i + 1) * 3]
        AA.append(codon_table[codon])
    return ''.join(AA)


def fasta_reader(filename):
    ''' read a fasta file. '''
    with open(filename, 'r') as file:
        name = file.readline().strip()
        seq = re.sub('[ \n\t]*', '', file.read())
        return (name, seq)


def parse_mutation(mutant_str, mutation_format):
    ''' parse a mutation description.
    e.g. 'c.123G>A' -> return (123,G,A)
     'p.R175H' -> return (175,R,H)
'''
    (regex, (ipos, ibefore, iafter)) = mutation_format
    result = regex.search(mutant_str)
    pos = int(result.groups()[ipos]) if result is not None else None
    before = result.groups()[ibefore] if result is not None else None
    after = result.groups()[iafter] if result is not None else None
    return (pos, before, after)


def convert(mutations_filename=os.path.join(__ROOT_DIR, 'mutation_data', 'L2L3database.csv'),
            seq_na_filename=os.path.join(__ROOT_DIR, 'seq_data', 'tp53_NM_000546.fasta'),
            seq_aa_filename=os.path.join(__ROOT_DIR, 'seq_data', 'P04637.fasta'),
            search_regions=list(regions.keys()),
            fout=sys.stdout,
            mode_use_mutation_column=False,
            filter_wildtype={},
            filter_only_for={},
            filter_omit={},
            columns_output_left=(),
            columns_output_right=(),
            columns_output_generated=(),
            columns_error=(),
            ):

    '''Converter.
    mutations_filename: mutation database (IARC TP53database, L2L3)
    seq_na_filename: coding sequence (only exons)
    seq_aa_filename: protein sequence
    search_regions: 変異の位置がここで指定した領域内である行のみが出力される
    fout: output
    ===DATABASE DEPEND PARAMETER===
    mode_use_mutation_column ==  True: column_mutation_NA,AAの変異情報を使う
                             == False: (codon_number,wildtype_codon,mutant_codon) の情報を使う
    filter_wildtype: wildtype行の指定．ここで指定したカラムが指定した文字列群である行は強制で出力
                     するとともに，配列も生成される
    filter_only_for: ここで指定したカラムが指定した文字列群である行だけが対象となる
    filter_omit: ここで指定したカラムが指定した文字列群である行は対象外とする
    columns_output_left: generatedの左側にそのまま出力する入力CSV上のカラム
    columns_output_right: generatedの右側にそのまま出力する入力CSV上のカラム
    columns_output_generated: 新しく生成するカラム
                              (領域名, mutated or wildtype, 変異塩基配列, 変異アミノ酸配列) の名前

    columns_error: エラー時に出力するカラム

    Usage:
    from genecoder import converter
    converter.convert(mutations_filename='***', *converter.T2T3)

    '''

    csv_writer = csv.writer(fout)

    # read the wildtype sequences
    (name_na, seq_na) = fasta_reader(seq_na_filename)
    (name_aa, seq_aa) = fasta_reader(seq_aa_filename)

    # output the header
    csv_writer.writerow(columns_output_left + columns_output_generated + columns_output_right)

    numerr_parse = 0
    numerr_wronginformation = 0
    numerr_wrongtranslation = 0
    with open(mutations_filename, 'rU') as mutations_file:
        try:
            for row in csv.DictReader(mutations_file):
                # L2L3database special (add mutation info)
                if mode_use_mutation_column is False and row[column_codon_number]:
                    L2L3_codonid = int(row[column_codon_number]) - 1
                    L2L3_wt = row[column_wildtype_codon]
                    L2L3_mutant = row[column_mutant_codon]
                    for w, m, p in zip(L2L3_wt, L2L3_mutant, range(3)):
                        if w is not m:
                            L2L3_description_na = (L2L3_codonid * 3 + p + 1, w, m)
                            L2L3_description_aa = (
                                codon_table[L2L3_wt], L2L3_codonid + 1, codon_table[L2L3_mutant])
                            break
                    row[column_mutation_NA] = 'c.{0}{1}>{2}'.format(
                        L2L3_description_na[0], L2L3_description_na[1], L2L3_description_na[2])
                    row[column_mutation_AA] = 'p.{0}{1}{2}'.format(
                        L2L3_description_aa[0], L2L3_description_aa[1], L2L3_description_aa[2])

                # filtering
                wildType = False
                for k, v in filter_wildtype.items():
                    if row[k] in v:
                        wildType = True
                        break
                if wildType is False:
                    ok = True
                    for k, v in filter_only_for.items():
                        if not row[k] in v:
                            ok = False
                            break
                    for k, v in filter_omit.items():
                        if row[k] in v:
                            ok = False
                            break
                    if not ok:
                        continue

                if wildType is False:
                    # get mutation information
                    try:
                        (na_pos, na_before, na_after) = parse_mutation(
                            row[column_mutation_NA], mutation_format_na)
                        (aa_pos, aa_before, aa_after) = parse_mutation(
                            row[column_mutation_AA], mutation_format_aa)
                        # convert an index format to python style (0 start)
                        na_pos -= 1
                        aa_pos -= 1
                    except Exception:
                        a_columns = ','.join(
                            ['{0}={1}'.format(k, v) for k, v in
                                zip(columns_error, (row[a_col] for a_col in columns_error))])
                        sys.stderr.write('error: [{0}] parse error. {1} or {2} is wrong.\n'.format(
                            a_columns, row[column_mutation_NA], row[column_mutation_AA]))
                        numerr_parse += 1
                        continue
                    # find a region including the mutation
                    for target_region_name in search_regions:
                        (left, right) = regions[target_region_name]
                        # convert an index format to python style
                        left -= 1
                        if left <= aa_pos < right:
                            # validate mutation
                            if na_before != seq_na[na_pos] or aa_before != seq_aa[aa_pos]:
                                a_columns = ','.join(
                                    ['{0}={1}'.format(k, v) for k, v in
                                        zip(columns_error,
                                            (row[a_col] for a_col in columns_error))])
                                sys.stderr.write(
                                    'error: [{0}] mutation information may be wrong.\n'.format(
                                        a_columns))
                                numerr_wronginformation += 1
                                break
                            # generate mutated sequence
                            seq_na_mutated = seq_na[:na_pos] + na_after + seq_na[na_pos + 1:]
                            seq_aa_mutated = seq_aa[:aa_pos] + aa_after + seq_aa[aa_pos + 1:]
                            seq_aa_mutated_transformed = NA_to_AA(
                                seq_na_mutated)[:-1]  # omit the last stop codon

                            # if seq_aa_mutated_transformed != seq_aa_mutated and not 'X' in
                            # seq_aa_mutated_transformed:
                            if seq_aa_mutated_transformed != seq_aa_mutated:
                                a_columns = ','.join(
                                    ['{0}={1}'.format(k, v) for k, v in
                                     zip(columns_error, (row[a_col] for a_col in columns_error))])
                                sys.stderr.write(
                                    'error: [{0}] translation may be wrong.\n'.format(a_columns))
                                numerr_wrongtranslation = 0
                                break
                            # --- output ---
                            (aa_left, aa_right) = regions[target_region_name]
                            # convert an index format to python style
                            aa_left -= 1
                            na_left, na_right = aa_left * 3, aa_right * 3
                            # left columns
                            output_items = []
                            output_items += [row[a_col] if row[a_col] !=
                                             '' else '' for a_col in columns_output_left]
                            # main columns
                            output_items += [target_region_name, 'mutated',
                                             seq_aa_mutated[aa_left:aa_right],
                                             seq_na_mutated[na_left:na_right]]
                            # right columns
                            output_items += [row[a_col] if row[a_col] !=
                                             '' else '' for a_col in columns_output_right]
                            # output the enter code
                            csv_writer.writerow(output_items)
                else:  # if wildType is True
                    for target_region_name in search_regions:
                        (aa_left, aa_right) = regions[target_region_name]
                        # convert an index format to python style
                        aa_left -= 1
                        (na_left, na_right) = (aa_left * 3, aa_right * 3)

                        # --- output ---
                        # left columns
                        output_items = []
                        output_items += [row[a_col] if row[a_col] !=
                                         '' else '' for a_col in columns_output_left]
                        # main columns
                        output_items += [target_region_name, 'wildtype',
                                         seq_aa[aa_left:aa_right], seq_na[na_left:na_right]]
                        # right columns
                        output_items += [row[a_col] if row[a_col] !=
                                         '' else '' for a_col in columns_output_right]
                        # output the enter code
                        csv_writer.writerow(output_items)

        except UnicodeDecodeError:
            sys.stderr.write('error')
    # output summary
    sys.stderr.write('errors: parse={0} wrong information={1} wrong translation={2}.\n'.format(
        numerr_parse, numerr_wronginformation, numerr_wrongtranslation))
