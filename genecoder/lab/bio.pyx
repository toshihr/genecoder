# cython: profile=True
from __future__ import absolute_import, division, print_function, unicode_literals
try:
    from future_builtins import map, zip
except ImportError:
    # Python 3 raise ImportError
    pass
import itertools
from fractions import Fraction

# codon table
codon_table = dict(zip((''.join(a_list) for a_list in itertools.product(
    'TCAG', repeat=3)), 'FFLLSSSSYYXXCCXWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))

# a nucleotide acid -> an element of GF(4)={0,1,alpha^1,alpha^2} table
NA2GF_table = dict()
GF2NA_table = dict()


def set_elements(s):
    global NA2GF_table
    global GF2NA_table
    NA2GF_table = dict(zip(s.upper() + s.lower(), [0, 1, 2, 3, 0, 1, 2, 3]))
    GF2NA_table = dict(zip(NA2GF_table.values(), NA2GF_table.keys()))

# global initialization
set_elements('ATGC')


def split_n(dat, n):
    ''' a generator to split the dat with length 1/n.

    example:
            >>> from genecoder.lab import bio
            >>> list(bio.split_n('atgcgcatg', 3))
            ['atg', 'cgc', 'atg']

            >>> list(bio.split_n([1,2,3,4,5,6], 3))
            [[1, 2, 3], [4, 5, 6]]

    '''
    assert len(dat) % 3 == 0

    for i in range(0, len(dat), n):
        yield dat[i:i + n]
    '''
    sio = io.StringIO(string)
    while True:
        s = sio.read(n)
        if s:
            yield s
        else:
            break
    '''


def get_codon12(s):
    ''' get a substring that contains only specific positions.

    example:
            >>> from genecoder.lab import bio
            >>> s = 'atgcgcatg'
            >>> str(bio.get_codon12(s))
            'atcgat'

    '''
    assert len(s) % 3 == 0
    return ''.join([a + b for a, b, in zip(s[0::3], s[1::3])])


def NA2GF(s, table=NA2GF_table):
    ''' convert a sequence to a polynomial.

    example:
            >>> from genecoder.lab import bio
            >>> s = 'atgcgcatg'
            >>> bio.NA2GF(s)
            [0, 1, 2, 3, 2, 3, 0, 1, 2]

    '''
    return list(map((lambda x: table[x]), s))


def GF2NA(p_x, capital=True, table=GF2NA_table):
    ''' convert a polynomial to a sequence.

    example:
            >>> from genecoder.lab import bio
            >>> a_x = [0, 1, 2, 3, 2, 3, 0, 1, 2]
            >>> str(bio.GF2NA(a_x))
            'ATGCGCATG'
            >>> str(bio.GF2NA(a_x, capital=False))
            'atgcgcatg'

    '''
    s = ''.join(map((lambda x: table[x]), p_x))
    return s.upper() if capital else s.lower()


def NA2AA(s):
    '''

            >>> from genecoder.lab import bio
            >>> str(bio.NA2AA('atgatg'))
            'MM'

    '''
    res = []
    for a_codon in split_n(s, 3):
        res.append(codon_table[a_codon.upper()])
    return str(''.join(res))


def codon_sort(s, n, k):
    ''' insert a parity strings to 3rd position.

    example:
            >>> from genecoder.lab import bio
            >>> s = 'AAAATCGGGGCT'
            >>> str(bio.codon_sort(s, n=6, k=4))
            'AATAACGGCGGT'

    '''
    assert n * 2 // 3 == k

    res = []
    for a_block in split_n(s, n):
        a_block_info = a_block[:k]
        a_block_parity = a_block[k:]
        a_block_shuffled = ''.join(
            [a + b + c for a, b, c in zip(a_block_info[0::2], a_block_info[1::2], a_block_parity)])
        res.append(a_block_shuffled)
    return str(''.join(res))


def get_similarity(s1, s2):
    '''

    '''
    assert len(s1) == len(s2)

    c = list(map((lambda x, y: x == y), s1, s2)).count(True)
    return Fraction(c, len(s1))


def RC(s1, s2):
    return get_similarity(NA2AA(s1), NA2AA(s2))
