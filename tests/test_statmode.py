# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals
from genecoder.main import main
import shlex
import os
import shutil
import csv

import pytest
pytest.main(__file__)

filename = 'L2L3seqs_event_T3T4.csv'
absname = os.path.join(os.path.dirname(__file__), 'data', filename)
outdir = os.path.join(os.path.dirname(__file__), '_tmpdir')
outfile = os.path.join(outdir, filename, 'STAT_ALL.csv')


class TestStatMode:

    @classmethod
    def setup_class(clazz):
        try:
            shutil.rmtree(outdir)
        except:
            pass

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

    def test_data(self):
        # run
        cmd = 'stat --coder n3_1 --gf4 atgc --outdir {0} --input {1}'.format(outdir, absname)
        assert main(shlex.split(cmd)) == 0

        with open(outfile, 'rU') as fin:
            reader = csv.DictReader(fin)
            line = next(reader)
            answer = {'region': 'L3', 'gf4': 'ATGC', 'n1': '104'}
            try:
                assert all([line[k] == v for k, v in answer.items()]) is True
            except AssertionError as e:
                print(line)
                raise e
