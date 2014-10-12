# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals
from genecoder.main import main
import shlex
import shutil
import os

filename = 'H7HARBD.csv'
absname = os.path.join(os.path.dirname(__file__), 'data', filename)
outdir = os.path.join(os.path.dirname(__file__), '_tmpdir')
outfile = os.path.join(outdir, 'tmp.fasta')


class TestAdd:

    @classmethod
    def setup_class(clazz):
        try:
            shutil.rmtree(outdir)
        except:
            pass
        os.mkdir(outdir)

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

    def test_simple1(self):
        cmd = 'csv2fasta 0 3 --input={0} --output={1}'.format(absname, outfile)
        assert main(shlex.split(cmd)) == 0

        with open(outfile, 'rU') as fin:
            header = fin.readline().strip()
            seq = fin.readline().strip()

            assert header[1:] == 'AAC40998 Influenza A virus (A/England/268/1996(H7N7))'
            assert seq == 'ggagcaaccagtgca'.upper()
