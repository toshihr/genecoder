from genecoder.main import main
import shlex

import pytest
pytest.main(__file__)


class TestDistance:

    @classmethod
    def setup_class(clazz):
        pass

    @classmethod
    def teardown_class(clazz):
        pass

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_simple1(self):
        cmd = 'distance --coder n3_1 --gf4 atgc --seq seq1:TTTCTTATTGTT'
        assert main(shlex.split(cmd)) == 0
