# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
try:
    from future_builtins import map
except ImportError:
    # Python 3 raise ImportError
    pass
import argparse
import sys
import csv
import configparser
import os
import itertools
from genecoder.lab import analyze
from genecoder.lab.fasta import Fasta
# the followings are necessary for parsing genecoder.ini
from genecoder.lab.codec import Coder_BCH, Coder_Convolution, coder_Iwadare, coder_SelfOrthogonal
from genecoder.lab.poly_gf4 import GF4
try:
    Coder_BCH
    Coder_Convolution
    coder_Iwadare
    coder_SelfOrthogonal
    GF4
except:
    pass

# --- global variables ---
parameters = {}
coder_list = {}


def loadResource(file_name):
    if not os.path.exists(file_name):
        raise Exception('{0} is not found.'.format(file_name))
    # A list of [information, n, k, coder]
    config = configparser.ConfigParser()
    config.read(file_name)
    section_name = 'coder'
    if section_name not in config.sections():
        raise Exception(('resource error. there is not the section "coder".'
                         'sections:{0}').format(config.sections()))
    a_coder_name_list = config.options(section_name)

    global coder_list
    coder_list = {}
    for a_coder_name in a_coder_name_list:
        coder_list[a_coder_name] = eval(config.get(section_name, a_coder_name))

    # coder_list = {
    # 	'bch1':('BCH code: g(x) = x + 1',3,2, Coder_BCH(g_x=[1,1],n=3)),
    # 	'bch2':('BCH code: g(x) = x^5 + a^2x^4 + a^2x^3 + a^2x^2 + a^2x + a',15,10,
    #       Coder_BCH(g_x=[1,3,3,3,3,2],n=15)),
    # 	'iwadare':('Iwadare code:',3,2, coder_Iwadare),
    # 	'selfortho':('self-orthogonal code:',3,2, coder_SelfOrthogonal),
    # }
    # print(coder_list)
    # sys.exit()


def set_parameters_CUI(argv):
    '''

    <yourscript> --file=outfile -q
    this code is same as the followings:
    <yourscript> -f outfile --quiet
    <yourscript> --quiet --file outfile
    <yourscript> -q -foutfile
    <yourscript> -qfoutfile
    '''

    # --- parse the options ---
    parser = argparse.ArgumentParser(description='This is a gene code program')
    # mode
    parser.add_argument('-showcoders', dest='mode', action='store_const',
                        const='showcoders', help='show the coders defined in genecoder.ini')
    parser.add_argument(
        '-stat', dest='mode', action='store_const', const='stat', help='do statistical tests')
    parser.add_argument('-stat_with_graph', dest='mode', action='store_const',
                        const='stat_with_graph', help='do statistical test and draw a graph')
    # now, solve mode is disabled
    # parser.add_argument('-solve'     , dest='mode', action='store_const', const='solve', type=int,
    # default=[], help='calculate a generator matrix from the sequences', nargs=1, metavar='n')
    # parameters
    parser.add_argument('-seq', dest='sequence', action='append', default=[],
                        help='add the DNA sequence to the analyze set', nargs=2,
                        metavar=('name', 'seq'))
    parser.add_argument('-code', dest='code', action='append',
                        default=[], help='add the code defined in the resource file')
    parser.add_argument('-gf4', dest='gf4', action='append', default=[
                        'ATGC'], help='add the coordinate of GF4 elements (0,1,2,3). e.g. ATGC')
    parser.add_argument('-gf4all', dest='gf4all', action='store_true',
                        default=False, help='use the all coordinates of GF4 elements')
    parser.add_argument('-all', dest='all', action='store_true',
                        default=False, help='use the all codes defined in the resource file')
    parser.add_argument('-fraction', dest='float', action='store_false',
                        default=True, help='the values are stored as a fraction style')
    parser.add_argument('-quiet', dest='verbose', action='store_false',
                        default=True, help='don\'t print status messages to stdout')
    parser.add_argument('-quieterr', dest='quieterr', action='store_true',
                        default=False, help='don\'t print error messages to stdout')
    parser.add_argument('-output', dest='output', nargs='?',
                        default='', help='output CSV file name or directory name in stat mode')
    parser.add_argument('-gui', dest='gui', action='store_true',
                        default=False, help='graphical user interface mode')
    parser.add_argument('infile', default='', nargs='?',
                        help=('FASTA file'
                              ' or CSV(comma separate) file in stat mode'), metavar='infile')
    options = parser.parse_args(argv)

    # --- set parameters ---
    global parameters

    parameters = {}

    def setParams_GF4():
        # set coordinate of GF4 elements
        parameters['GF4'] = []
        if options.gf4 != []:
            parameters['GF4'] = list(
                set(parameters['GF4'] + list(map(lambda x: x.upper(), options.gf4))))
        if options.gf4all:
            parameters['GF4'] = [''.join(x)
                                 for x in itertools.permutations('ATGC')]

    def setParams_coders():
        # set coders
        if options.all is False:
            parameters['coders'] = list(
                map((lambda x: x.lower()), options.code))
        else:
            parameters['coders'] = [x for x in coder_list.keys()]

    def checkParams_coders():
        # coder check
        for coder_exec in parameters['coders']:
            if coder_exec not in coder_list.keys():
                print('code {0} is not defined!'.format(coder_exec))
                sys.exit(1)
        if not parameters['gui mode'] and len(parameters['coders']) == 0:
            print('choose more than one coder!')
            sys.exit(1)

    # initialize options for gui mode
    parameters['input database'] = os.path.abspath(options.infile) if options.infile != '' else ''
    parameters['output stat'] = os.path.abspath(options.output) if options.output != '' else ''
    parameters['output'] = os.path.abspath(options.output) if options.output != '' else ''
    parameters['input sequences'] = []

    # set other options
    analyze.fracStyle = options.float
    analyze.verboseWarning = options.verbose
    analyze.quietErr = options.quieterr
    parameters['gui mode'] = options.gui
    parameters['mode'] = options.mode if options.mode is not None else 'standard'

    if parameters['mode'] == 'stat_with_graph':
        parameters['mode'] = 'stat'
        parameters['drawGraph'] = True
    else:
        parameters['drawGraph'] = False

    if parameters['mode'] == 'showcoders':
        # show coders mode
        pass
    elif parameters['mode'] == 'stat':
        # set coders
        setParams_coders()
        # set coordinate of GF4 elements
        setParams_GF4()
        # --- check ---
        if len(parameters['output stat']) == 0:
            print('output is needed in stat mode.')
            sys.exit(1)
        if len(parameters['input database']) == 0:
            print('infile is needed in stat mode.')
            sys.exit(1)
        # coder check
        checkParams_coders()

    else:
        # standard mode
        # set sequences
        if options.infile != '':
            # read a fasta file
            fasta = Fasta(open(options.infile))
            parameters['input sequences'] += list(fasta.blocks)
        if options.sequence != []:
            # append the tuple (name, seq)
            parameters['input sequences'] += options.sequence
        # set coders
        setParams_coders()
        # set coordinate of GF4 elements
        setParams_GF4()
        # --- check ---
        # input sequence check
        if not parameters['gui mode'] and len(parameters['input sequences']) == 0:
            print('Wrong usage. try -h option to describe the usage.')
            sys.exit(1)
        # coder check
        checkParams_coders()


def analyze_and_output():
    # --- statistical tests mode ---
    if parameters['mode'] == 'stat':
        analyze.analyze_survivalTest_for_database(
            out_dir=parameters['output stat'],
            database=parameters['input database'],
            parameters=parameters,
            coder_list=coder_list)
    else:
        try:
            if parameters['output'] != '':
                output = open(parameters['output'], 'w', newline='')
            else:
                output = sys.stdout

            csv_writer = csv.writer(output, delimiter=',', quotechar='"')

            # --- calculate generator matrix mode ---
            if parameters['mode'] == 'solve':
                analyze.analyze_estimate_code(csv_writer, parameters)
            # --- standard analyze mode ---
            elif parameters['mode'] == 'standard':
                # --- output the header ---
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
                    'coder_detail')
                csv_writer.writerow(header)
                # --- analyze ---
                for a_coordinate_of_GF4 in parameters['GF4']:
                    for (coder_id, coder_detail, n, k, a_coder) in map(
                            (lambda x: [x, ] + list(coder_list[x])), parameters['coders']):
                        for (a_name, s1, s2, AA1, AA2, RC, simirarity) in analyze.gen_RC_distance(
                                seqs=parameters['input sequences'],
                                coder=a_coder,
                                GF4_coordinate=a_coordinate_of_GF4):
                            a_line = [
                                a_name, s1, s2, AA1, AA2, RC, simirarity, a_coordinate_of_GF4,
                                coder_id, coder_detail]
                            csv_writer.writerow(a_line)
        except IOError as e:
            raise e
        else:
            if parameters['output'] != '':
                output.close()


def run(argv):
    # argv = sys.argv[1:] + ['-solve','3','-seq','seq1','TTTCTTATTGTT']
    # argv = sys.argv[1:] +
    #        ['-solve','9','-solve','6','-solve','3','-seq','seq1',s1,'-seq','seq2','TTTCTTATTGTT']
    # argv = sys.argv[1:] + ['-solve','9','-seq','seq1',s1]
    set_parameters_CUI(argv)

    if parameters['mode'] == 'showcoders':
        print('--- support coders defined in genecoder.ini ---')
        print(list(coder_list.keys()))
        print()
    elif parameters['mode'] == 'stat':
        analyze_and_output()
    elif parameters['mode'] == 'solve':
        analyze_and_output()
    else:
        if parameters['gui mode']:
            from PySide import QtGui
            from genecoder.gui import MyMainWindow
            app = QtGui.QApplication(sys.argv[:1])
            dlg = MyMainWindow()
            dlg.ok = False
            dlg.show()
            if app.exec_() != 0:
                return 1
            if dlg.ok is False:
                return 0
        else:
            analyze_and_output()

    return 0


# --- global initialization ---
loadResource(os.path.join(os.path.dirname(__file__), 'genecoder.ini'))


def main():
    if len(sys.argv) == 1:
        return run(['-h'])
    else:
        return run(sys.argv[1:])


if __name__ == '__main__':
    sys.exit(main())
