#!/usr/bin/python3
import pstats
import cProfile
import shlex
# Note that:
# tottime: この関数に費した総時間で、この関数から呼び出された関数は含まない
# cumtime: この関数に費した総時間で、この関数から呼び出された関数を含む
#
# References:
# [1] http://docscythonja.zouri.jp/src/tutorial/profiling_tutorial.html

# cmdline = ('-output {0} -quiet -quieterr -all -gf4all -seq aa atgcatgcatgc'
#            ' -seq bb cgatatgcatgcatgc').format(os.devnull)
cmdline = ('-output aaa -quiet -quieterr -stat -code bch_n3_1 -code bch_n3_2'
           ' test/database_L2L3_test.csv')
args = shlex.split(cmdline)
cProfile.runctx('genecoder.run(args)', globals(), locals(), 'Profile.prof')

s = pstats.Stats('Profile.prof')
s.strip_dirs().sort_stats('time').print_stats()
