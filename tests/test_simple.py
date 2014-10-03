from genecoder.main import run
import shlex


class TestAdd:

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
        assert run(shlex.split('-code bch_n3_1 -seq seq1 TTTCTTATTGTT')) == 0
