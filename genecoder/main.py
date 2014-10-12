# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals
import sys
import csv
from genecoder.resource import CODERS
from genecoder.lab import analyze
from genecoder import cli


def mode_distance(args):
    header = (
        'name',
        'dna sequence',
        'enocoded dna sequence',
        'protein',
        'encoded protein',
        'RC',
        'similarity between original dna sequence and encoded dna sequence',
        'coordinate of GF(4) elements',
        'coder_id',
        'coder_detail'
    )

    if args['--output'] is None:
        fout = sys.stdout
    else:
        fout = open(args['--output'], 'w')
    csv_writer = csv.writer(fout)

    csv_writer.writerow(header)

    for a_coordinate_of_GF4 in args['--gf4']:
        for (coder_id, coder_detail, n, k, a_coder) in map(
                (lambda x: [x, ] + list(CODERS[x])), args['--coder']):
            for (a_name, s1, s2, AA1, AA2, RC, simirarity) in analyze.gen_RC_distance(
                    seqs=args['input sequences'],
                    coder=a_coder,
                    GF4_coordinate=a_coordinate_of_GF4):
                a_line = [
                    a_name, s1, s2, AA1, AA2, RC, simirarity, a_coordinate_of_GF4,
                    coder_id, coder_detail]
                csv_writer.writerow(a_line)

    return 0


def mode_stat(args):
    analyze.analyze_survivalTest_for_database(out_dir=args['--outdir'], database=args['--input'],
                                              target_gf4=args['--gf4'],
                                              target_coders=args['--coder'],
                                              drawGraph=args['--graph'])
    return 0


def mode_list(args):
    from genecoder.resource import CODERS
    print('--- support coders defined in genecoder.ini ---')
    print(list(CODERS.keys()))
    print()
    return 0


def mode_gui(args):
    from PySide import QtGui
    from genecoder import gui

    gui.pipe_args = args
    app = QtGui.QApplication(sys.argv[:1])
    dlg = gui.MyMainWindow()
    dlg.ok = False
    dlg.show()
    if app.exec_() != 0:
        return 1
    if dlg.ok is False:
        return 0


def mode_csv2fasta(args):
    from genecoder.csv2fasta import csv2fasta

    if args['--input'] is None:
        fin = sys.stdin
    else:
        fin = open(args['--input'], 'rU')

    if args['--output'] is None:
        fout = sys.stdout
    else:
        fout = open(args['--output'], 'w')

    csv2fasta(fin, fout, args['<idx_name>'], args['<idx_seq>'], args['<length>'], to_upper=True)

    return 0


def main(argv=sys.argv[1:]):

    args = cli.main(argv)

    analyze.fracStyle = args['--fraction']
    analyze.verboseWarning = args['--quiet-verbose'] or args['--quiet']
    analyze.quietErr = args['--quiet']

    if args['list']:
        return mode_list(args)
    elif args['stat']:
        return mode_stat(args)
    elif args['distance']:
        return mode_distance(args)
    elif args['gui']:
        return mode_gui(args)
    elif args['csv2fasta']:
        return mode_csv2fasta(args)
    return 0


if __name__ == '__main__':
    sys.exit(main())
