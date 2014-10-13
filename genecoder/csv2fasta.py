# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals
import csv
import textwrap


def csv2fasta(fin, fout, idx_name, idx_seq, read_length=None, to_upper=True):
    '''read csv, then return multiple fasta.

        if read_length is None then read full length

    '''

    csv_reader = csv.reader(fin)

    # skip csv header
    next(csv_reader)

    for a_line in csv_reader:
        if all(a_symbol in 'atgcATGC' for a_symbol in a_line[idx_seq]):
            length = len(a_line[idx_seq]) if read_length is None else read_length
            # skip null fasta line
            if length == 0:
                continue

            seq = textwrap.fill(a_line[idx_seq][:length], 80) + '\n'
            # output
            fout.write('>' + a_line[idx_name] + '\n')
            if to_upper:
                fout.write(seq.upper())
            else:
                fout.write(seq)
