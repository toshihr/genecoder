# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals
from genecoder.main import run
import shlex
import os
import shutil
import csv

filename = 'L2L3seqs_event_T3T4.csv'
absname = os.path.join(os.path.dirname(__file__), 'data', filename)
outdir = os.path.join(os.path.dirname(__file__), '_tmpdir')
outfile = os.path.join(outdir, filename, 'STAT_ALL.csv')


class TestStatMode:

    @classmethod
    def setup_class(clazz):
        try:
            shutil.rmtree(outdir)
        except FileNotFoundError:
            pass

    @classmethod
    def teardown_class(clazz):
        try:
            shutil.rmtree(outdir)
        except FileNotFoundError:
            pass

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_data(self):
        # run
        assert run(shlex.split('-stat -code bch_n3_1 -o {0} {1}'.format(outdir, absname))) == 0

        with open(outfile, 'r', newline='') as fin:
            reader = csv.DictReader(fin, delimiter=',', quotechar='"')
            line = next(reader)
            answer = {'region': 'L3', 'gf4': 'ATGC', 'n1': '104'}
            assert all([line[k] == v for k, v in answer.items()]) is True
