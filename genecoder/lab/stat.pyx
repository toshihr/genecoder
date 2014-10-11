# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals
import math
from scipy.stats import norm
from scipy.stats import chi2
import numpy as np
import matplotlib as mpl
mpl.use('Qt4Agg')
import matplotlib.pyplot as plt
import collections

'''
生存分析について
線形回帰：連続量
ロジスティック回帰：（０、１）など２値しかとらない変数
Cox回帰（比較ハザードモデル）：観察期間と観察終了時のイベント発生の有無

中途打ち切り例（censored cases）

1. Kaplan-Meier method (product limit method)
少数例でも正確に生存率を計算できる。ケースが多いときにも用いられるが、その場合はCutler-Ederer methodが使用できる。


2. generalized wilcoxon
References:
http://d.hatena.ne.jp/ryamada22/20120331/1333009284
http://aoki2.si.gunma-u.ac.jp/R/Gen-Wil.html
http://aoki2.si.gunma-u.ac.jp/lecture/Survival/gw.html

# group: 0 は A 群，1 は B 群を表す (青木先生のサイトでは1がA群, 2がB群)
# event: 1 は 死亡，0 は 生存

data:  time: time, event: event, group: group
W = 63.0000, Var(W) = 1320.3126, Z = 1.7338, p-value = 0.08295


3. logrank test, generalized wilcoxon test
logrank test: 生存数が少なくなる観察期間の後半の観察結果に左右される傾向が強い
generalized wilcoxon test: センサリングのパターンに依存する度合いが強い
'''


def plot_km(data_list, legend=['legend1', 'legend2', 'legend3'], color=True, label=(
        'time', 'cumulative probability'), title=('Kaplan-Meier', 24), pvalues=[], outfile=''):
    # 日本語表示
    # import matplotlib.font_manager as fm
    # prop = fm.FontProperties(fname='/usr/share/fonts/truetype/fonts-japanese-gothic.ttf')
    # plt.xlabel('あ', fontsize=20, fontproperties=prop)

    if color:
        # axis1にWildTypeが含まれているときを想定. axis1の一部がMutationTypeとなる
        styles = ['Dk-', 'or-', 'sg-']
        styles_tag = ['or', 'or', 'sg']
    else:
        styles = ['Dk-', 'ok--', 'sk-.']
        styles_tag = ['ok', 'ok', 'sk']

    plt.clf()

    for a_data, a_style, a_style_tag, a_label in zip(data_list, styles, styles_tag, legend):
        if not a_data:
            continue
        pre_x = 0
        pre_y = 1.0
        # レジェンド描画対策のため最初の一点を描画する
        plt.plot([pre_x, pre_x], [pre_y, pre_y], a_style, label=a_label)

        for x, y, tag in zip(a_data[0], a_data[3], a_data[5] if len(a_data) >= 6 else np.array(
                [False] * a_data[0].size)):
            if pre_y == y:
                # 平行移動のとき
                # 横線
                plt.plot([pre_x, x], [pre_y, y], a_style[1:])
            else:
                # 右下に落ちるとき
                # 横線
                plt.plot([pre_x, x], [pre_y, pre_y], a_style[1:])
                # 縦線
                plt.plot([x, x], [pre_y, y], a_style[1:])
            # 点
            plt.plot([x, x], [y, y], a_style[:2] if not tag else a_style_tag)
            pre_x = x
            pre_y = y

    plt.title(title[0], size=title[1])
    plt.ylim(0.0, 1.1)
    plt.xlabel(label[0], fontsize=20)
    plt.ylabel(label[1], fontsize=20)

    # legend
    # plt.legend()

    # pvalue
    for a_pvalue, index in zip(pvalues, range(len(pvalues))):
        plt.text(0.65, 0.95 - 0.05 * index, a_pvalue, fontsize=12, transform=plt.gca().transAxes)

    if outfile == '':
        plt.show()
    else:
        plt.savefig(outfile)


def km(time_list, event_list, tag_list=None):
    ''' Kaplan-Meier method
    input:
    time_list: 生存時間
    event_list: イベント(死亡)発生なら1, それ以外(打ち切り)は0

    return:
    a tuple (numpy.array time, flag, probability, cumulative probability, SE)

    References:
    http://aoki2.si.gunma-u.ac.jp/R/km-surv.html
    http://aoki2.si.gunma-u.ac.jp/lecture/Survival/k-m.html

    >>> from genecoder.lab.stat import *
    >>> g = np.array([0,0,1,1,0,0,1,1,0,1,1,0,0,0,1,1,1,0,0,1,0,1,0,1,0,1,1,0,0,0])
    >>> e = np.array([1,0,1,1,0,0,1,0,0,0,1,1,0,0,0,1,0,0,1,1,0,1,0,0,1,0,1,0,0,0])
    >>> t = np.array([2,20,5,1,3,17,2,3,15,14,12,13,11,11,10,8,8,3,7,3,6,2,5,4,2,3,1,3,2,1])
    >>> g_a = g==0
    >>> g_b = g!=0
    >>> data1 = km(t[g_a],e[g_a])
    >>> data2 = km(t[g_b],e[g_b])
    >>> for a_line in zip(*data1):
    ...     print( ','.join(list(map(lambda x: '{0:>6.3}'.format(x) if isinstance(x, float) else
    ...  '{0:>6}'.format(x), a_line))) )
         1,     1,   1.0,   1.0,   nan
         2,     0, 0.933, 0.933,0.0622
         2,     0, 0.929, 0.867,0.0847
         2,     1,   1.0, 0.867,   nan
         3,     1,   1.0, 0.867,   nan
         3,     1,   1.0, 0.867,   nan
         3,     1,   1.0, 0.867,   nan
         5,     1,   1.0, 0.867,   nan
         6,     1,   1.0, 0.867,   nan
         7,     0, 0.857, 0.743, 0.129
        11,     1,   1.0, 0.743,   nan
        11,     1,   1.0, 0.743,   nan
        13,     0,  0.75, 0.557, 0.169
        15,     1,   1.0, 0.557,   nan
        17,     1,   1.0, 0.557,   nan
        20,     1,   1.0, 0.557,   nan

    各カラムはtime,trunc,p,P,SEを表す.
    trunc=1(STOP)->Pは変化しない
    trunc=0(EVENT)->Pが変化する(グラフの高さが変わる)
    同じ時刻でかつtrunc=1(STOP)の点は上書きされるのでみえない

    #>>> plot_km((data1,data2))

    '''
    assert len(time_list) == len(event_list)

    if tag_list is None:
        data = np.array([time_list, event_list])
    else:
        assert len(tag_list) == len(event_list)
        data = np.array([time_list, event_list, tag_list])

    # 欠損値対策
    data = data[:, ~np.isnan(data).any(0)]

    n = len(data[0])

    # timeでソート(イベントと打ち切りの生存期間が同じ場合の細工として、打ち切りの方が後になるよう打ち切りデータの生存時間に微小値を加える)
    if np.sum(data[0] > 0) > 0:
        fraction = np.amin(data[0][data[0] > 0]) / 1000
        modified_time = np.where(data[1], data[0], data[0] + fraction)

        order = np.argsort(modified_time)
        data[0] = data[0][order]
        data[1] = data[1][order]
        if tag_list is not None:
            data[2] = data[2][order]

    # 整数ベクトル
    I = np.arange(1, n + 1)
    # 生存確率 (打ち切り(event = 0)なら1)
    p = np.where(data[1] == 0, 1, (n - I) / (n - I + 1))
    # 累積確率
    P = np.cumprod(p)
    # 標準誤差 (イベント時点でのみ計算される)
    se = data[1] / (n - I + 1) / (n - I + 1)
    SE = np.where(data[1], P * np.sqrt(np.cumsum(se)), np.nan)

    if tag_list is None:
        result = (data[0], 1 - data[1], p, P, SE)
    else:
        result = (data[0], 1 - data[1], p, P, SE, data[2])

    return result


def logrank(group_list, event_list, time_list, method='SAS'):
    '''

    method: SAS or Tominaga

    return:
            a tuple of X-squared,p-value,method='SAS' or 'Tominaga',table

    >>> from genecoder.lab.stat import *
    >>> g = [0,0,1,1,0,0,1,1,0,1,1,0,0,0,1,1,1,0,0,1,0,1,0,1,0,1,1,0,0,0]
    >>> e = [1,0,1,1,0,0,1,0,0,0,1,1,0,0,0,1,0,0,1,1,0,1,0,0,1,0,1,0,0,0]
    >>> t = [2,20,5,1,3,17,2,3,15,14,12,13,11,11,10,8,8,3,7,3,6,2,5,4,2,3,1,3,2,1]
    >>> result = logrank(g,e,t)
    >>> print('chi={0:.5},p-value={1:.5},method={2:.5}'.format(result[0],result[1],result[2]))
    chi=3.3805,p-value=0.065971,method=SAS

    #chi=3.38053228738,p-value=0.0659707459488,method=SAS
    >>> print(np.array_str(result[3],max_line_width=255,precision=4))
    [[  0.       2.       2.      16.      14.      30.       0.0667   1.0667   0.9333]
     [  2.       2.       4.      15.      12.      27.       0.1481   2.2222   1.7778]
     [  0.       1.       1.      12.      10.      22.       0.0455   0.5455   0.4545]
     [  0.       0.       0.       9.       7.      16.       0.       0.       0.    ]
     [  0.       1.       1.       9.       6.      15.       0.0667   0.6      0.4   ]
     [  0.       0.       0.       8.       5.      13.       0.       0.       0.    ]
     [  1.       0.       1.       7.       5.      12.       0.0833   0.5833   0.4167]
     [  0.       1.       1.       6.       5.      11.       0.0909   0.5455   0.4545]
     [  0.       0.       0.       6.       3.       9.       0.       0.       0.    ]
     [  0.       0.       0.       6.       2.       8.       0.       0.       0.    ]
     [  0.       1.       1.       4.       2.       6.       0.1667   0.6667   0.3333]
     [  1.       0.       1.       4.       1.       5.       0.2      0.8      0.2   ]
     [  0.       0.       0.       3.       1.       4.       0.       0.       0.    ]
     [  0.       0.       0.       3.       0.       3.       0.       0.       0.    ]
     [  0.       0.       0.       2.       0.       2.       0.       0.       0.    ]
     [  0.       0.       0.       1.       0.       1.       0.       0.       0.    ]]

    reference:
    http://aoki2.si.gunma-u.ac.jp/R/logrank.html
    '''

    assert len(group_list) == len(event_list) == len(time_list)

    # 欠損値対策
    data = np.array([group_list, event_list, time_list], dtype=float)
    data = data[:, ~np.isnan(data).any(0)]

    return logrank_numpy(data=data, method=method)


def logrank_numpy(data, method='SAS'):
    ''' numpyの行列としてデータを受け取りlogrank検定を行う関数.
    data = numpy.array([group, event, time], dtype=float)
    欠損値はあらかじめ省いておく必要がある
    破壊的関数(ソートするため)
    '''
    # ソート(time昇順)
    order = np.argsort(data[2])
    data[0] = data[0][order]
    data[1] = data[1][order]
    data[2] = data[2][order]

    g_set = np.unique(data[0])
    e_set = np.unique(data[1])
    t_set = np.unique(data[2])
    assert len(g_set) == 2 and len(e_set) == 2

    # (g,e,t) -> freq
    # freq = collections.Counter([tuple(data[:, i]) for i in range(0, len(data[0]))])
    # t -> freq
    # note that the length of each keys() is not the same each other
    # one should be access by the elements of t_set
    freq_g1e0 = collections.Counter(data[2, (data[0] == g_set[0]) * (data[1] == e_set[0])])
    freq_g1e1 = collections.Counter(data[2, (data[0] == g_set[0]) * (data[1] == e_set[1])])
    freq_g2e0 = collections.Counter(data[2, (data[0] == g_set[1]) * (data[1] == e_set[0])])
    freq_g2e1 = collections.Counter(data[2, (data[0] == g_set[1]) * (data[1] == e_set[1])])
    tg_g1e0 = np.array([freq_g1e0[t] for t in sorted(t_set)])  # tg[,1]
    tg_g1e1 = np.array([freq_g1e1[t] for t in sorted(t_set)])  # tg[,2]
    tg_g2e0 = np.array([freq_g2e0[t] for t in sorted(t_set)])  # tg[,3]
    tg_g2e1 = np.array([freq_g2e1[t] for t in sorted(t_set)])  # tg[,4]

    k = len(t_set)
    nia = np.sum(tg_g1e0) + np.sum(tg_g1e1)
    nib = np.sum(tg_g2e0) + np.sum(tg_g2e1)
    na = np.zeros(k)
    na[0] = nia
    na[1:] = (np.ones(k) * nia - np.cumsum(tg_g1e0 + tg_g1e1))[:-1]
    nb = np.zeros(k)
    nb[0] = nib
    nb[1:] = (np.ones(k) * nib - np.cumsum(tg_g2e0 + tg_g2e1))[:-1]
    da = tg_g1e1
    db = tg_g2e1
    dt = da + db
    nt = na + nb
    d = dt / nt
    O = np.array([np.sum(da), np.sum(db)])
    ea = na * d
    eb = nb * d
    E = np.array([np.sum(ea), np.sum(eb)])

    result = np.column_stack((da, db, dt, na, nb, nt, d, ea, eb))
    # print(na)
    # result = np.column_stack((da,db,dt))
    # result = (da,db,dt,na,nb,nt,d,ea,eb)

    old_settings = np.seterr(all='ignore')
    if method == 'SAS':
        v = np.nansum(dt * (nt - dt) / (nt - 1) * na / nt * (1 - na / nt))
        temp = (np.sum(da) - np.sum(na * d))
        chi = temp * temp / v
    else:
        temp = O - E
        chi = np.sum(temp * temp / E)
    P = 1.0 - chi2.cdf(chi, 1)
    np.seterr(**old_settings)

    return (chi, P, method, result)


def gen_wilcox(group_list, event_list, time_list):
    '''

    return:
            a tuple of W,var(W),z-value,p-value

    >>> from genecoder.lab.stat import *
    >>> g = [0,0,1,1,0,0,1,1,0,1,1,0,0,0,1,1,1,0,0,1,0,1,0,1,0,1,1,0,0,0]
    >>> e = [1,0,1,1,0,0,1,0,0,0,1,1,0,0,0,1,0,0,1,1,0,1,0,0,1,0,1,0,0,0]
    >>> t = [2,20,5,1,3,17,2,3,15,14,12,13,11,11,10,8,8,3,7,3,6,2,5,4,2,3,1,3,2,1]
    >>> gen_wilcox(g,e,t)
    (63, 1320.3126436781608, 1.7338126143696353, 0.082951336727237868)

    '''

    assert len(group_list) == len(event_list) == len(time_list)

    # 欠損値対策
    data = np.array([group_list, event_list, time_list], dtype=float)
    data = data[:, ~np.isnan(data).any(0)]

    return gen_wilcox_numpy(data)


def gen_wilcox_numpy(data):
    ''' logrankと同様numpy行列に対して一般化wilcox検定を行う. '''
    def get_U(tx, sx, ty, sy):
        if (tx < ty and sx == 1 and sy == 1) or (tx <= ty and sx == 1 and sy == 0):
            return -1
        elif (tx > ty and sx == 1 and sy == 1) or (tx >= ty and sx == 0 and sy == 1):
            return 1
        else:
            return 0

    # ソート(group昇順)
    order = np.argsort(data[0])
    data[0] = data[0][order]
    data[1] = data[1][order]
    data[2] = data[2][order]

    g_set = np.unique(data[0])
    e_set = np.unique(data[1])
    assert len(g_set) == 2 and len(e_set) == 2

    n = len(data[0])
    na = sum(data[0] == g_set[0])
    nb = n - na

    W = 0
    for i in range(0, na):
        for j in range(na, n):
            W += get_U(data[2, i], data[1, i], data[2, j], data[1, j])

    W_var = 0
    for i in range(0, n):
        temp = 0
        for j in range(0, n):
            temp += get_U(data[2, i], data[1, i], data[2, j], data[1, j])
        W_var += temp * temp
    W_var = na * nb * W_var / n / (n - 1)
    try:
        Z = W / math.sqrt(W_var)
        P = (1.0 - norm.cdf(abs(Z))) * 2
    except ZeroDivisionError:
        Z = np.nan
        P = np.nan
    return (W, W_var, Z, P)
