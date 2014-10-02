#!/usr/bin/python
import csv
import sys

if len(sys.argv) < 3:
    print(('Usage: csv2fasta [idx of name] [idx of seq] [length] < [input csv]'
           ' > [output fasta]. initial idx 0.'))
    quit()

idx_name = int(sys.argv[1])
idx_seq = int(sys.argv[2])
if len(sys.argv) > 3:
    dna_length = int(sys.argv[3])
else:
    dna_length = 0

csv_reader = csv.reader(sys.stdin, delimiter=',', quotechar='"')
# read a header
header = next(csv_reader)
# treat each line
for a_line in csv_reader:
    # check
    if all(a_symbol in 'atgcATGC' for a_symbol in a_line[idx_seq]):
        if dna_length == 0:
            dna_length = len(a_line[idx_seq])
        print(">" + a_line[idx_name])
        print(a_line[idx_seq][:dna_length].upper())
