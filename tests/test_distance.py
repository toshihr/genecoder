from genecoder.main import main
import shlex
import os
import shutil

import pytest
pytest.main(__file__)

outdir = os.path.join(os.path.dirname(__file__), '_tmpdir')
outfile = os.path.join(outdir, 'tmp.csv.gz')


class TestDistance:

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
        cmd = 'distance --coder n3_1 --gf4 atgc --seq seq1:TTTCTTATTGTT'
        assert main(shlex.split(cmd)) == 0

    def test_output_csv_gzip(self):
        cmd = 'distance --coder n3_1 --gf4 atgc --seq seq1:TTTCTTATTGTT --output ' + outfile
        assert main(shlex.split(cmd)) == 0
        assert os.path.isfile(outfile) is True
