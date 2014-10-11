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

        >>> from genecoder.lab.poly_gf4 import *
        >>> p = [1,0,1]
        >>> str(Poly.to_string(p))
        '1x^2+0x^1+1'

        >>> p = [1,2,3]
        >>> str(Poly.to_string(p))
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

        >>> from genecoder.lab.poly_gf4 import *
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

        >>> from genecoder.lab.poly_gf4 import *
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
