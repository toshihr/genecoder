# -*- coding: utf-8 -*-

'''genecoder

usage:
    genecoder distance [--quiet-verbose | --quiet] [--fraction]
                       [--coder=<coder> ... | --all] [--gf4=<gf4> | --gf4all]
                       [--output <output>] (--seq=<label_na> ... | --input=<fasta>)
    genecoder stat     [-quiet | -Quiet] [--fraction]
                       [--coder=<coder> | --all] [--gf4=<gf4> | --gf4all]
                       [--graph] (--outdir <outdir>) (--input <dataset>)
    genecoder list
    genecoder gui
    genecoder -h | --help
    genecoder -v | --version

options:
    -q, --quiet-verbose             don't print verbose messages
    -Q, --quiet                     don't print verbose & error messages
    -f, --fraction                  store the values as a fraction style
    -c=<coder>, --coder=<coder>     add coder
    -a, --all                       add all supported codes
    --gf4=<gf4>                     set the correspondings with the elements of GF(4)
    --gf4all                        use all combinations of correspondings
    --seq=<label_na>                add 'label:nucleotide sequence' style sequence
    -i=<fasta>, --input=<fasta>     add sequences from the FASTA file
    -o=<output>, --output=<output>  result CSV file
    --outdir=<outdir>               result directory
    -h, --help                      show this help
    -v, --version                   show version

'''

from __future__ import absolute_import, division, print_function
try:
    from future_builtins import filter
except ImportError:
    # Python 3 raise ImportError
    pass
from docopt import docopt
import sys
import re
from genecoder.resource import NAME, VERSION, CODERS
from genecoder.lab.fasta import Fasta


def tune_args(args):
    # coders
    if args['--all']:
        args['--coder'] = list(CODERS.keys())

    # input sequences (distance mode)
    if args['distance']:
        input_sequences = []
        # --seq: convert ['label:na',] to [('label', 'na'),]
        for label_na in args['--seq']:
            try:
                label, na = label_na.split(':')
                input_sequences.append((label, na))
            except ValueError:
                raise ValueError('--seq {0} is wrong usage'.format(args['--seq']))
        # --input: read FASTA file
        if len(args['--input']) > 0:
            fasta = Fasta(open(args['--input'][0], 'rU'))
            input_sequences.extend(list(fasta.blocks))
        # add new entry
        args['input sequences'] = input_sequences

    return args


def validate_args(args):
    # validate coder
    set_selected = set(args['--coder'])
    set_all = set(CODERS.keys())
    if set_selected.issubset(set_all) is False:
        raise ValueError('coder name is wrong: {0}'.format(set_selected.difference(set_all)))
    # validate input sequences
    seqs = [seq for label, seq in args['input sequences']]
    p = re.compile('^[atgcATGC]+$')
    not_passed = list(filter(lambda s: p.match(s) is None, seqs))
    if len(not_passed) > 0:
        raise ValueError('sequence contains unknown symbol: {0}'.format(not_passed))


def main(argv=sys.argv[1:]):
    # parse argv. no options, with -v, -h will execute sys.exit()
    args = docopt(__doc__, argv=argv, version='{0} {1}'.format(NAME, VERSION), options_first=False)
    # tune
    args = tune_args(args)
    # validate
    validate_args(args)

    print(args)


if __name__ == '__main__':
    sys.exit(main())
