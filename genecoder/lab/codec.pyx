# -*- coding: utf-8 -*-
# cython: profile=True
from __future__ import absolute_import, division, print_function, unicode_literals
import abc
from genecoder.lab.poly_gf4 import GF4, Poly


def poly2str(poly, base='x'):
    cdef int l
    cdef int i

    l = len(poly) - 1
    s = []
    for i in range(len(poly)):
        if poly[i] > 1:
            coeff = '{0}'.format(poly[i])
        elif poly[i] == 1:
            coeff = ''
        else:
            continue

        if l - i > 1:
            s.append('{0}{1}^{2}'.format(coeff, base, l - i))
        elif l - i == 1:
            s.append('{0}{1}'.format(coeff, base))
        else:
            # const
            if poly[i] > 0:
                s.append('{0}'.format(poly[i]))

    if len(s) == 0:
        s.append('0')

    return ' + '.join(s)


class Coder(object, metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def encode(self, a_x):
        pass
#    @abc.abstractmethod
#    def decode(self, a_x): pass


class Coder_Cyclic(Coder):

    ''' Cyclic code (including BCH code, Reed-Solomon code) over the field.
    parameter:
            g_x: a generator polynomial with the degree n-k (should be a primitive polynomial)
            n: a length of code
            field: a finite field such as GF2, GF4
            k: a length of information
    algorithm:
            ref. [1] p.154
    '''

    def __init__(self, g_x, n, field=GF4):
        self.g_x = g_x
        self.n = n
        self.field = field
        self.k = n - (len(g_x) - 1)

    def encode_1block(self, a_x):
        ''' Cyclic code over the field.
        input:
                a_x is a data polynomial to encode with the degree k
        example:
                (n=7,k=4) Cyclic code over GF(4)
                >>> from genecoder.lab.codec import *
                >>> g_x = [1,0,1,1] # x^3 + x + 1
                >>> dat = [1,1,0,1]
                >>> coder = Coder_Cyclic(g_x, n=7)
                >>> coder.encode_1block(dat)
                [1, 1, 0, 1, 0, 0, 1]

        '''
        assert len(a_x) == self.k, 'len(a_x)={0},self.k={1}'.format(len(a_x), self.k)

        b_x = Poly.shift(a_x, self.n - self.k)
        r_x = Poly.mod(b_x, self.g_x, self.field)
        return Poly.plus(b_x, r_x)

    def encode(self, a_x):
        '''Cyclic code over the field.
        input:
                a_x is a data polynomial to encode with the degree k
        example:
                (n=7,k=4) Cyclic code over GF(4)
                >>> from genecoder.lab.codec import *
                >>> from genecoder.lab.poly_gf4 import GF2
                >>> g_x = [1,0,1,1] # x^3 + x + 1
                >>> dat = [1,1,0,1,1,1,0,1]
                >>> coder = Coder_Cyclic(g_x, n=7, field=GF2)
                >>> coder.encode(dat)
                [1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1]

        '''
        res = []
        cdef int i
        for i in range(0, len(a_x), self.k):
            res.extend(self.encode_1block(a_x[i:i + self.k]))
        return res

    def __str__(self):
        return 'n = {0}, k = {1}, g(x) = {2}'.format(self.n, self.k, poly2str(self.g_x))


class Coder_Convolution(Coder):

    ''' convolutional code over the field.
    parameter:
            G: a generator matrix (transfer function matrix) with the dimension (k,n)
                            should be stored in column major e.g. [G11,G21,G12,22,G13,G23]
            n: a length of code
            field: a finite field such as GF2, GF4
            k: a length of information
    algorithm:
            ref. [2] p.250
    example:
            2/3 convolutional code over GF(4)
            >>> from genecoder.lab.codec import *
            >>> G13 = [1,0,1] # D^2 + 1
            >>> G23 = [1,1] # D + 1
            >>> G = [[1,],[0,],[0,],[1,],G13,G23]
            >>> coder = Coder_Convolution(G=G, n=3, k=2, field=GF4)
            >>> coder.encode([0,1,2,3,0,1,2,3,0,1])
            [0, 1, 1, 2, 3, 0, 0, 1, 2, 2, 3, 2, 0, 1, 2]

    '''

    def __init__(self, G, n, k, field=GF4):
        self.G = G
        self.n = n
        self.k = k
        self.field = field
        assert len(G) == k * n

    def encode(self, a_x):
        def M_i_t(i, t):
            ''' return the M_i(t) i.e. the ith column of information block with time t.'''
            return a_x[self.k * t + i] if 0 <= t else 0

        # block loop
        W = []
        for t in range(0, len(a_x) // self.k):
            W_t = []
            # calculate W_j(t) j=0,1,...,n-1
            for j in range(0, self.n):
                # W_j(t) = sum(M_i(t,t-1,...) * G_ij ; i=0,1,...,k-1)
                # 出力のブロックtにおける第j列目の値は、入力のブロック群G_ijの第i列目の値の和をi=0,...,k-1まで足し合わせたものとなる
                W_j_t = 0
                for i in range(0, self.k):
                    G_ij = self.G[j * self.k + i]
                    for s in range(0, len(G_ij)):
                        if G_ij[-1 - s] != 0:
                            W_j_t = self.field.plus(W_j_t, M_i_t(i, t - s))
                W_t.append(W_j_t)
            W.extend(W_t)

        return W

    def __str__(self):
        return 'n = {0}, k = {1}, G = [{2}]'.format(
            self.n, self.k, ','.join(map(lambda p: poly2str(p, 'D'), self.G)))


class Coder_Linear(Coder):

    ''' linear code over the field.
    parameter:
            G: a generator matrix (transfer function matrix) with the dimension (k,n)
                            should be stored in column major e.g. [G11,G21,G12,22,G13,G23]
            n: a length of code
            field: a finite field such as GF2, GF4
            k: a length of information
    algorithm:
    example:
            (3,2)-linear code over GF(4)
            >>> from genecoder.lab.codec import *
            >>> G = [1,0,0,1,1,1]
            >>> coder = Coder_Linear(G=G, n=3, k=2, field=GF4)
            >>> coder.encode([0,1,2,3,0,1,2,3,0,1])
            [0, 1, 1, 2, 3, 1, 0, 1, 1, 2, 3, 1, 0, 1, 1]

    '''

    def __init__(self, G, n, k, field=GF4):
        self.G = G
        self.n = n
        self.k = k
        self.field = field
        assert len(G) == k * n

    def encode(self, a_x):
        W = []
        for t in range(0, len(a_x) // self.k):
            W_t = []
            # calculate W_j(t) j=0,1,...,n-1
            for j in range(0, self.n):
                # W_j(t) = sum(M_i(t,t-1,...) * G_ij ; i=0,1,...,k-1)
                W_j_t = 0
                for i in range(0, self.k):
                    G_ij = self.G[j * self.k + i]
                    a_x_i = a_x[t * self.k + i]
                    a_x_i_mul_G_ij = self.field.product(a_x_i, G_ij)
                    W_j_t = self.field.plus(W_j_t, a_x_i_mul_G_ij)
                W_t.append(W_j_t)
            W.extend(W_t)

        return W

    def __str__(self):
        return 'G = [{0}]'.format(','.join(self.G))

# === famous code ===
# Iwadare code (n=3, field=GF(4))
# Ref: [2] p.272 オリジナルの岩垂符号は2元に限定
# G13 = D^7 + D^5
# G23 = D^4 + D^3
coder_Iwadare = Coder_Convolution(
    G=[[1, ], [0, ], [0, ], [1, ], [1, 0, 1, 0, 0, 0, 0, 0], [1, 1, 0, 0, 0]], n=3, k=2, field=GF4)
# self-orthogonal code (r=2/3, field=GF(4)) same as Wyner-Ash code
# Ref: [2] p.277
# G13 = D^2 + 1
# G23 = D+1
coder_SelfOrthogonal = Coder_Convolution(
    G=[[1, ], [0, ], [0, ], [1, ], [1, 0, 1], [1, 1]], n=3, k=2, field=GF4)
