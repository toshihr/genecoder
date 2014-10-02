# -*- coding: utf-8 -*-
# cython: profile=True
from __future__ import absolute_import, division, print_function, unicode_literals
try:
    from future_builtins import map, zip
except ImportError:
    # Python 3 raise ImportError
    pass
import abc
import itertools
from genecoder.lab.genericmatrix import GenericMatrix

''' library of the polynomials over GF(2^2).
:すべてのメソッドは非破壊的
:Note that the polynomial's coefficients are in decreasing powers.
教科書における多項式表現とは逆順であることに注意。高位次数の情報部が左側、下位次数のパリティ部が右側にくる。
すなわち[x^n-1,x^n-2,...,x^2,x^1,x^0]=[information;parity]
左右反転すれば教科書における多項式表現と一致する。

Note:
GF(q)は、qを法とする整数の演算で作ることはできない。そこでGF(p)上の多項式を用いる。
P(x)をGF(p)上のm次の既約多項式とする。
P(x)は素体GF(p)における素数pに似た役割を演じる。
GF(p^m)が構成できる。
GF(p^m)は色々な作り方ができるができあがるものは本質的には同じもの。元の順番をかえれば同じ演算表となる。

GF(p^m)上の演算について
[2]p.56
乗除：べき表現を用いると、
a^i * a^j = a^{(i+j)mod (p^m-1)}
a^i / a^j = a^{(i-j)mod (p^m-1)}
と簡単に行える。
加減：多項式基底による展開、またはベクトル表現に変換して行う必要がある。
注意 GF(p^m)上の多項式の係数をmod (p^m)するのは危険。元の整数への対応順によって結果が変わるため。あくまで演算表にしたがうべき。

GF(4) = GF(2^2)の演算表:
(GF(2)上の2次の既約多項式はx^2+x+1しかない)
+		0:0		1:1		2:x		3:x+1
0:0		0:0		1:1		2:x		3:x+1
1:1		1:1		0:0		3:x+1	2:x
2:x		2:x		3:x+1	0:0		1:1
3:x+1	3:x+1	2:x		1:1		0:0

*		0:0		1:1		2:x		3:x+1
0:0		0:0		0:0		0:0		0:0
1:1		0:0		1:1		2:x		3:x+1
2:x		0:0		2:x		3:x+1	1:1
3:x+1	0:0		3:x+1	1:1		2:x

・べき表現との対応: 0=0 1=1 x=a x+1=a^2
・本ライブラリ上の表記との対応: 0=0=(00) 1=1=(01) a=2=(10) a^2=3=(11)
・2進2次元ベクトルだと思えば、和差はbitwise XORで計算できる. ATGCの対応づけは別の話なので問題ない.
  次数nの多項式の係数ベクトルはP={0,1}^2n次元ベクトルで表現できる. つまり長さ2nの2進配列
  あるいは、長さnの整数配列を用い、演算をXORにすれば整数を2進表現すれば係数ベクトルになっているため同じこと.
  積除演算はこううまくはいかない
GF(2^2) = Z_2(x)/<x^2+x+1>
Reed-Solomon：リードソロモン符号はGF(2^p)

References:
[1] 藤原良, 神保雅一, 符号と暗号の数理, 共立出版
[2] 今井秀樹, 符号理論
'''
# === finite fields ===


class GF(object, metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def plus(f1, f2):
        pass

    @abc.abstractmethod
    def minus(f1, f2):
        pass

    @abc.abstractmethod
    def product(f1, f2):
        pass

    @abc.abstractmethod
    def divide(f1, f2):
        pass

    @abc.abstractmethod
    def to_string(f):
        pass


class GF2(GF):

    @staticmethod
    def plus(f1, f2):
        return (f1 + f2) % 2

    @staticmethod
    def minus(f1, f2):
        return (f1 - f2) % 2

    @staticmethod
    def product(f1, f2):
        return f1 * f2

    @staticmethod
    def divide(f1, f2):
        return f1 // f2

    @staticmethod
    def to_string(f):
        return str(f)


class GF4(GF):

    '''
    GF(2^2) = {0,1,x,x+1}
    0->0, 1->1, 2->x, 3->x+1
    '''
    # static variables
    plus_table = dict(zip(itertools.product((0, 1, 2, 3), repeat=2),
                          (0, 1, 2, 3,
                           1, 0, 3, 2,
                           2, 3, 0, 1,
                           3, 2, 1, 0)))
    product_table = dict(zip(itertools.product((0, 1, 2, 3), repeat=2),
                             (0, 0, 0, 0,
                              0, 1, 2, 3,
                              0, 2, 3, 1,
                              0, 3, 1, 2)))
    inverse_table = {
        1: 1,
        2: 3,
        3: 2,
    }
    to_string_table = {0: '0', 1: '1', 2: 'a', 3: 'a^2'}
    minus_table = {}
    divide_table = {}

    def __init__(self):
        pass

    @staticmethod
    def init_static_variables():
        for (a, b), c in GF4.plus_table.items():
            GF4.minus_table[(c, b)] = a
        for (a, b), c in GF4.product_table.items():
            if b != 0:
                GF4.divide_table[(c, b)] = a

    @staticmethod
    def plus(f1, f2):
        return GF4.plus_table[(f1, f2)]

    @staticmethod
    def minus(f1, f2):
        return GF4.minus_table[(f1, f2)]

    @staticmethod
    def product(f1, f2):
        return GF4.product_table[(f1, f2)]

    @staticmethod
    def divide(f1, f2):
        return GF4.divide_table[(f1, f2)]

    @staticmethod
    def to_string(f):
        return GF4.to_string_table[f]

GF4.init_static_variables()

# === polynomials ===


class Poly:

    @staticmethod
    def fix(a_x):
        ''' omit the first terms which have zero coefficient. '''
        a_x = a_x[:]
        while a_x[0] == 0 and len(a_x) > 1:
            a_x.pop(0)
        return a_x

    @staticmethod
    def product(a_x, b, field=GF4):
        ''' a(x) * b
        a(x): a polynomial over GF(4)
        b: an element of GF(4)
        '''
        return [field.product(f, b) for f in a_x]

    @staticmethod
    def shift(a_x, n):
        a_x = a_x[:]
        ''' x^nを掛ける '''
        assert n >= 0
        a_x += [0, ] * n
        return a_x

    @staticmethod
    def unshift(a_x, n):
        ''' x^nで割る '''
        assert n >= 0
        return a_x[0:-n]

    @staticmethod
    def to_string(a_x, field=GF4):
        ''' return the usual representation.

        >>> from lab.poly_gf4 import *
        >>> p = [1,0,1]
        >>> Poly.to_string(p)
        '1x^2+0x^1+1'

        >>> p = [1,2,3]
        >>> Poly.to_string(p)
        '1x^2+ax^1+a^2'

        '''
        terms = []
        d = len(a_x) - 1
        if d > 0:
            for f in a_x[0:-1]:
                terms.append('{0}x^{1}'.format(field.to_string(f), d))
                d -= 1
        terms.append(field.to_string(a_x[-1]))
        return '+'.join(terms)

    @staticmethod
    def mod(a_x, b_x, field=GF4):
        ''' a(x) mod b(x) over the field.

        >>> from lab.poly_gf4 import *
        >>> a = [3,0,0] # a^2x^2
        >>> b = [1,0,1] # x^2 + 1
        >>> Poly.mod(a,b,field=GF4)
        [3]

        '''
        assert field == GF2 or field == GF4
        assert b_x[0] > 0
        a_x = a_x[:]
        if len(a_x) < len(b_x):
            return a_x
        else:
            m = field.divide(a_x[0], b_x[0])
            mb_x = Poly.product(b_x, m)
            # r(x) = a(x) - m * b(x)
            r_x = []
            for i in range(len(mb_x)):
                t = field.minus(a_x.pop(0), mb_x[i])
                if not (t == 0 and len(r_x) == 0):
                    r_x.append(t)
            # add a constant term
            r_x.extend(a_x)
            return Poly.mod(r_x, b_x)

    @staticmethod
    def plus(a_x, b_x, field=GF4):
        ''' a(x) + b(x) over the field

        >>> from lab.poly_gf4 import *
        >>> a_x =   [1,0,1,1] #       x^3 + x + 1
        >>> b_x = [1,1,0,1,1] # x^4 + x^3 + x + 1
        >>> Poly.plus(a_x, b_x, field=GF4)
        [1, 0, 0, 0, 0]

        '''
        # d_min = min(len(a_x), len(b_x))
        d_max = max(len(a_x), len(b_x))
        new_list = []
        for i in range(0, d_max):
            a = a_x[-1 - i] if i < len(a_x) else 0
            b = b_x[-1 - i] if i < len(b_x) else 0
            new_list.insert(0, field.plus(a, b))
        return new_list


class Matrix:

    ''' a matrix over the finite field.
    this class use the genericmatrix.py at http://www.mit.edu/~emin/source_code/py_ecc/

    >>> from lab.poly_gf4 import *
    >>> m = Matrix(3,3,[0,3,2]+[1,0,3]+[2,1,0],field=GF4)
    >>> m
    <matrix
     0 3 2
     1 0 3
     2 1 0>
    >>> m+m
    <matrix
     0 0 0
     0 0 0
     0 0 0>
    >>> m*m
    <matrix
     0 2 2
     1 0 2
     1 1 0>
    >>> m.toLine()
    '032 103 210'

    '''

    def __init__(self, rows, cols, mat=None, field=GF4):
        self.rows = rows
        self.cols = cols
        self.field = GF4
        self.data = GenericMatrix((self.rows, self.cols), zeroElement=0, identityElement=1,
                                  add=field.plus, mul=field.product, sub=field.minus,
                                  div=field.divide)
        if type(mat) == GenericMatrix:
            for r in range(0, self.rows):
                self.data.SetRow(r, mat.GetRow(r))
        elif mat is not None:
            self.set(mat)

    def __repr__(self):
        return repr(self.data)

    def __setitem__(self, i_j, data):
        '''__setitem__((i,j),data) sets item row i and column j to data.'''
        (i, j) = i_j
        self.data[i, j] = data

    def __getitem__(self, i_j):
        '''__getitem__((i,j)) gets item at row i and column j.'''
        (i, j) = i_j
        return self.data[i, j]

    def __mul__(self, other):
        return Matrix(self.rows, other.cols, self.data * other.data, self.field)

    def __add__(self, other):
        return Matrix(self.rows, self.cols, self.data + other.data, self.field)

    def __sub__(self, other):
        return Matrix(self.rows, self.cols, self.data - other.data, self.field)

    def __hash__(self):
        return hash(tuple(flatten(self.data.data)))

    def __eq__(self, other):
        return self.data.data == other.data.data

    def inverse(self):
        return Matrix(self.rows, self.cols, self.data.Inverse(), self.field)

    def set(self, list_row_major):
        ''' set [a11,a12,...,a1m,a21,...,a2m,...,anm] '''
        for r in range(0, self.rows):
            self.data.SetRow(r, list_row_major[self.cols * r:self.cols * (r + 1)])

    def toLine(self):
        rows = [self.data.GetRow(i) for i in range(0, self.rows)]
        return ' '.join(list(map(lambda r: ''.join([str(x) for x in r]), rows)))

    def solve(self, b):
        ''' Solve the equation Ax = b, return a column vector x.
        b: a list (a column vector)
        '''
        assert self.rows == len(b)
        M = self.data.MakeSimilarMatrix((self.rows, self.cols + 1), 'z')
        for r in range(0, self.rows):
            M.SetRow(r, self.data.GetRow(r) + [b[r], ])

        num_rank = 0
        for c in range(0, self.cols):
            # 列cの中で対角成分から検索をはじめ、0でない行を見つける
            r_id = M.FindRowLeader(num_rank, c)
            if r_id >= 0:
                # その成分を1 にする
                M.MulRow(r_id, GF4.divide(1, M[r_id, c]))
                # 階段の角になる行と交換する
                if num_rank != r_id:
                    M.SwapRows(num_rank, r_id)
                # 階段の角となる行以外を0にする
                for r in itertools.chain(range(0, num_rank), range(num_rank + 1, self.rows)):
                    # M[r,] += -M[r,c] * M[num_rank,]
                    M.MulAddRow(GF4.minus(0, M[r, c]), num_rank, r)
                num_rank += 1
        if num_rank == self.cols:
            return M.GetColumn(self.cols)[:num_rank]
        elif num_rank < self.cols:
            return 'INFINITY rank={0},cols={1}'.format(num_rank, self.cols)
        else:
            raise ValueError('can not solved')

    def rank(self):
        '''
        >>> from lab.poly_gf4 import *
        >>> M = Matrix(3,3,[0,3,2]+[1,0,3]+[2,1,0],field=GF4)
        >>> M.rank()
        3

        >>> M = Matrix(3,3,[2,1,0]+[0,0,2]+[0,0,1],field=GF4)
        >>> M.rank()
        2

        '''
        M = self.data.Copy()
        num_rank = 0
        for c in range(0, self.cols):
            # 列cの中で対角成分から検索をはじめ、0でない行を見つける
            r_id = M.FindRowLeader(num_rank, c)
            if r_id >= 0:
                # その成分を1 にする
                M.MulRow(r_id, GF4.divide(1, M[r_id, c]))
                # 階段の角になる行と交換する
                if num_rank != r_id:
                    M.SwapRows(num_rank, r_id)
                # 階段の角となる行以外を0にする
                for r in itertools.chain(range(0, num_rank), range(num_rank + 1, self.rows)):
                    # M[r,] += -M[r,c] * M[num_rank,]
                    M.MulAddRow(GF4.minus(0, M[r, c]), num_rank, r)
                num_rank += 1
        return num_rank


def is_iter(l):
    try:
        iter(l)
        return True
    except TypeError:
        return False


def flatten(xs):
    if '__iter__' in dir(xs):
        for ys in map(flatten, iter(xs)):
            for y in ys:
                yield y
    else:
        yield xs


def calc_generator_from_codes(x_list, u_list, n, k):
    ''' calculate the generator G from the equation X=UG.
    x_list: a list of codes with length n
    u_list: a list of codes with length k
    the sizes of these list should be k

    X = (x_1,...,x_n) where x_i = (x_list[0][i],...,x_list[k-1][i]), i=1,...,n is column vectors
    U = (u_1,...,u_k) where u_i = (u_list[0][i],...,u_list[k-1][i]), i=1,...,k is column vectors
    G = (g_1,...,g_n)
    X=UG is same as the following:
    x_1 = Ug_1
    ...
    x_n = Ug_n

    >>> from lab.poly_gf4 import *
    >>> x_list = [[1,1,1],[0,1,1]]
    >>> u_list = [[1,1]  ,[0,1]  ]
    >>> calc_generator_from_codes(x_list, u_list, 3, 2)
    <matrix
     1 0 0
     0 1 1>

    >>> x_list = [[1,1,1],[1,1,1]]
    >>> u_list = [[1,1]  ,[1,1]  ]
    >>> calc_generator_from_codes(x_list, u_list, 3, 2)
    Traceback (most recent call last):
         ...
    ValueError: codes may not be independent

    '''
    assert len(x_list) == len(u_list)
    assert len(x_list) == k
    # X: (k,n)-matrix
    X = Matrix(k, n, mat=list(flatten(list(x_list))), field=GF4)
    # U: (k,k)-matrix
    U = Matrix(k, k, mat=list(flatten(list(u_list))), field=GF4)
    # G = X * U.inverse
    try:
        U.inverse()
    except ValueError:
        raise ValueError('codes may not be independent')

    return U.inverse() * X


def _test():
    import doctest
    doctest.testmod(verbose=False)

if __name__ == '__main__':
    _test()
