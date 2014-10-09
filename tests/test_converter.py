# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals
from genecoder.converter import convert, L2L3, IARC
import os
import shutil
import csv

database_L2L3 = os.path.join(os.path.dirname(__file__),
                             '../genecoder/rawdata/mutation_data/L2L3database.csv')
database_IARC = os.path.join(os.path.dirname(__file__),
                             '../genecoder/rawdata/IARC-TP53DATABASE/TP53SomaticR17.csv')

outdir = os.path.join(os.path.dirname(__file__), '_tmpdir')
outfile = os.path.join(outdir, 'tmp.csv')


class TestStatMode:

    @classmethod
    def setup_class(clazz):
        try:
            shutil.rmtree(outdir)
        except:
            pass
        os.makedirs(outdir)

    @classmethod
    def teardown_class(clazz):
        try:
            shutil.rmtree(outdir)
        except:
            pass

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_L2L3(self):
        with open(outfile, 'w') as fout:
            convert(mutations_filename=database_L2L3, fout=fout,
                    search_regions=['L2', 'L3'], **L2L3)

        with open(outfile, 'rU') as fin:
            reader = csv.DictReader(fin)
            line = next(reader)
            answer = {'PubMed': '11289122', 'OS(months)': '29.4'}
            try:
                assert all([line[k] == v for k, v in answer.items()]) is True
            except AssertionError as e:
                print(line)
                raise e

    def test_IARC(self):
        with open(outfile, 'w') as fout:
            convert(mutations_filename=database_IARC, fout=fout, **IARC)

        with open(outfile, 'rU') as fin:
            reader = csv.DictReader(fin)
            line = next(reader)
            answer = {'Mutation_ID': '1561', 'Short_topo': 'COLON'}
            try:
                assert all([line[k] == v for k, v in answer.items()]) is True
            except AssertionError as e:
                print(line)
                raise e
