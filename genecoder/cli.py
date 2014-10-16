# -*- coding: utf-8 -*-

'''code analysis for genes

usage:
    genecoder distance [--quiet-verbose | --quiet] [--fraction] [--compress ]
                       (--coder=<coder> ... | --all)
                       (--gf4=<gf4> ... | --gf4all)
                       [--output=<output>]
                       (--seq=<label_na> ... | --input=<fasta>)
    genecoder stat     [-quiet | -Quiet] [--fraction]
                       (--coder=<coder> | --all)
                       (--gf4=<gf4> | --gf4all)
                       [--graph] (--outdir=<outdir>) (--input=<dataset>)
    genecoder gui      [--fraction]
                       [--coder=<coder> ... | --all]
                       [--gf4=<gf4> ... | --gf4all]
                       [--graph] [--output=<output> | --outdir=<outdir>]
                       [--input=<fasta> | --input=<dataset>]
    genecoder list
    genecoder csv2fasta <idx_name> <idx_seq> [<length>] [--input=<csv>]
                        [--output=<output>]
    genecoder -h | --help
    genecoder -v | --version

modes:
    distance    calculate RC distance
    stat        survival tests
    list        show supported coders

options:
    -q, --quiet-verbose             don't print verbose messages
    -Q, --quiet                     don't print verbose & error messages
    -f, --fraction                  store the values as a fraction style
    -c=<coder>, --coder=<coder>     add coder
    -a, --all                       add all supported codes
    -g, --graph                     draw graph
    --compress                      output zipped csv, works with --output
    --gf4=<gf4>                     element correspondings, e.g. ATGC, ACGT
    --gf4all                        use all combinations of correspondings
    --seq=<label_na>                add 'label:nucleotide sequence'
    -i=<fasta>, --input=<fasta>     add sequences from the FASTA file
    -o=<output>, --output=<output>  result CSV file
    --outdir=<outdir>               result directory
    -h, --help                      show this help
    -v, --version                   show version

references:
- Sato Keiko, Toshihide Hara, and Masanori Ohya. "The code structure of the p53
  DNA-binding domain and the prognosis of breast cancer patients." Bioinformati
  cs 29.22 (2013): 2822-2825.

'''

from __future__ import absolute_import, division, print_function, unicode_literals
try:
    from future_builtins import filter
except ImportError:
    # Python 3 raise ImportError
    pass
from docopt import docopt
import sys
import re
import itertools
from genecoder.resource import CODERS
from genecoder.lab.fasta import Fasta
from genecoder.version import VERSION


def tune_args(args):
    '''
        - if '--add' is True then set --coder list to use all coders.
        - add 'input sequences' entry based on '--seq' & '--input' options.
    '''
    # input
    if len(args['--input']) > 0:
        args['--input'] = args['--input'][0]
    else:
        args['--input'] = None

    # coders: '--all'
    if args['--all']:
        args['--coder'] = list(CODERS.keys())
        args['--all'] = False

    # GF4: '--gf4', '--gf4all'
    args['--gf4'] = list(map(lambda s: s.upper(), args['--gf4']))
    if args['--gf4all']:
        args['--gf4'].extend([''.join(x) for x in itertools.permutations('ATGC')])
        args['--gf4all'] = False
    args['--gf4'] = list(set(args['--gf4']))

    # input sequences (distance mode): '--seq', '--input'
    args['input sequences'] = []
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
        if args['--input']:
            try:
                fasta = Fasta(open(args['--input'], 'rU'))
            except IOError as e:
                print(e)
                sys.exit(1)
            input_sequences.extend(list(fasta.blocks))
        # add new entry
        args['input sequences'] = input_sequences

    # csv2fasta mode
    args['<idx_name>'] = int(args['<idx_name>']) if args['<idx_name>'] else None
    args['<idx_seq>'] = int(args['<idx_seq>']) if args['<idx_seq>'] else None
    args['<length>'] = int(args['<length>']) if args['<length>'] else None

    return args


def validate_args(args):
    p = re.compile('^[atgcATGC]+$')
    # validate coder
    set_selected = set(args['--coder'])
    set_all = set(CODERS.keys())
    if set_selected.issubset(set_all) is False:
        raise ValueError('coder name is wrong: {0}'.format(set_selected.difference(set_all)))
    # validate input sequences
    seqs = [seq for label, seq in args['input sequences']]
    not_passed = list(filter(lambda s: p.match(s) is None, seqs))
    if len(not_passed) > 0:
        raise ValueError('sequence contains unknown symbol: {0}'.format(not_passed))
    # validate GF4
    not_passed = list(filter(lambda s: len(s) != 4 or p.match(s) is None, args['--gf4']))
    if len(not_passed) > 0:
        raise ValueError('GF4 is wrong: {0}'.format(not_passed))


def main(argv=sys.argv[1:]):
    # parse argv. no options, with -v, -h will execute sys.exit()
    args = docopt(__doc__, argv=argv, version='genecoder {0}'.format(VERSION), options_first=False)

    if __name__ == '__main__':
        print(args)

    # tune
    args = tune_args(args)
    # validate
    validate_args(args)

    if __name__ == '__main__':
        print(args)

    return args


if __name__ == '__main__':
    sys.exit(main())
