# -*- coding: utf-8 -*-
# cython: profile=True
from __future__ import absolute_import, division, print_function, unicode_literals
try:
    from future_builtins import map, zip
except ImportError:
    # Python 3 raise ImportError
    pass
import os
import shutil
import csv
import tarfile
import numpy as np
from collections import OrderedDict
from itertools import compress
from genecoder.lab import bio
from genecoder.lab import stat
from genecoder.resource import CODERS

# global parameters
FracStyle = True
verboseWarning = True
quietErr = False


def pack(absdir, erase=False):
    # .tar.gz -> mode=gz
    arcName = absdir + '.tar.bz2'
    th = tarfile.open(arcName, 'w:bz2')
    for root, dirs, files in os.walk(absdir):
        for f in files:
            fullpath = os.path.join(root, f)
            relativepath = fullpath[
                len(os.path.dirname(absdir) + os.path.sep):]
            # archive relative path
            th.add(fullpath, arcname=relativepath)
    th.close()

    if erase:
        shutil.rmtree(absdir, ignore_errors=True)

    return arcName


def gen_table(csv_reader, columns=None, accept_filters=None, omit_filters={}):
    ''' 指定した列だけをcsv_readerから読み込み、行のリストを返すジェネレータ関数.
    出力されるカラムの順番はcolumnsの順(CSVファイル上の順ではない！)
    columns: ここで指定された列の情報が読み込まれる. Noneの場合すべての列が読み込まれる
    accept_filters: これが指定された場合、すべての指定されたキーが指定した値だったときのみ読み込む
    omit_filters: どれかひとつでも指定されたキーが指定した値だったときは読み込まない
    omit_filtersがaccept_filtersより優先される

    example:
    >>> import csv
    >>> from genecoder.lab import analyze
    >>> fromDB = 'id,c1,c2,c3\\n001,1,2,3\\n002,2,3,4'
    >>> list(analyze.gen_table(csv.reader(fromDB.split('\\n'))))
    [['001', '1', '2', '3'], ['002', '2', '3', '4']]
    >>> list(analyze.gen_table(csv.reader(fromDB.split('\\n')), columns=['id','c1']))
    [['001', '1'], ['002', '2']]
    >>> list(analyze.gen_table(csv.reader(fromDB.split('\\n')),
    ... columns=['id','c1'], accept_filters={'c2':'3'}))
    [['002', '2']]
    >>> list(analyze.gen_table(csv.reader(fromDB.split('\\n')),
    ... columns=['id','c1'], omit_filters={'c2':'3'}))
    [['001', '1']]

    '''
    # read a header
    header = next(csv_reader)
    if columns is None:
        columns = header
    # make an index list
    index_columns = [header.index(a_column) for a_column in columns]
    index_value__accept_filters = dict(
        [(header.index(k), v) for k, v in accept_filters.items()]) if accept_filters else None
    index_value__omit_filters = dict(
        [(header.index(k), v) for k, v in omit_filters.items()])
    # treat each line
    for a_line in csv_reader:
        # omit filters
        skip = any(
            a_line[k] == v for k, v in index_value__omit_filters.items())
        if skip:
            continue
        # accept filters
        skip = not all(a_line[k] == v for k, v in index_value__accept_filters.items(
        )) if accept_filters else False
        if skip:
            continue
        # make a trimmed line, then append to the results
        yield [a_line[a_index] for a_index in index_columns]


def gen_columns_from_table(table, column_indexes,
                           index_value__accept_filters=None, index_value__omit_filters={}):
    ''' テーブルの指定した列を返すジェネレータ関数.
    accept_filters: {index:value}
    omit_filters: {index:value}
    example:
    >>> import csv
    >>> from genecoder.lab import analyze
    >>> fromDB = [['001', '1', '2', '3'], ['002', '2', '3', '4']]
    >>> list(analyze.gen_columns_from_table(fromDB, [0,1]))
    [['001', '1'], ['002', '2']]
    >>> list(analyze.gen_columns_from_table(fromDB, 0))
    [['001'], ['002']]
    >>> list(analyze.gen_columns_from_table(fromDB, [0,1], {1:'2'}))
    [['002', '2']]
    >>> list(analyze.gen_columns_from_table(fromDB, [0,1], None, {1:'2'}))
    [['001', '1']]

    '''
    if not (type(column_indexes) == list or type(column_indexes) == tuple):
        column_indexes = [column_indexes]
    flags_columns = [k in column_indexes for k in range(len(table[0]))] if len(
        table) > 0 else []
    for a_line in table:
        # omit filters
        skip = any(
            a_line[k] == v for k, v in index_value__omit_filters.items())
        if skip:
            continue
        # accept filters
        skip = not all(a_line[k] == v for k, v in index_value__accept_filters.items(
        )) if index_value__accept_filters else False
        if skip:
            continue
        # make a trimmed line, then append to the results
        yield list(compress(a_line, flags_columns))


def survivalTest(data, km=True):
    ''' 生存分析を行う.
    input:
     data: numpy.array([group,event,time],dtype=float)
    output: 結果(logrank, generalized wilcox, Kaplan-Meier)を返す. 検定ができない場合はNoneを返す
    破壊的関数
    example:
    >>> import numpy as np
    >>> from genecoder.lab import analyze
    >>> g = [0,0,1,1,0,0,1,1,0,1,1,0,0,0,1,1,1,0,0,1,0,1,0,1,0,1,1,0,0,0]
    >>> e = [1,0,1,1,0,0,1,0,0,0,1,1,0,0,0,1,0,0,1,1,0,1,0,0,1,0,1,0,0,0]
    >>> t = [2,20,5,1,3,17,2,3,15,14,12,13,11,11,10,8,8,3,7,3,6,2,5,4,2,3,1,3,2,1]
    >>> stat_wilcox,stat_logrank,stat_km_g1,stat_km_g2 = analyze.survivalTest(np.array([g,e,t],
    ... dtype=float))
    >>> print(stat_km_g1[3]) # 累積確率
    [ 1.          0.93333333  0.86666667  0.86666667  0.86666667  0.86666667
      0.86666667  0.86666667  0.86666667  0.74285714  0.74285714  0.74285714
      0.55714286  0.55714286  0.55714286  0.55714286]
    >>> print(stat_km_g2[3]) # 累積確率
    [ 0.92857143  0.85714286  0.78571429  0.71428571  0.64285714  0.64285714
      0.64285714  0.64285714  0.53571429  0.42857143  0.42857143  0.42857143
      0.21428571  0.21428571]
    >>> print(np.array_str(stat_logrank[3],max_line_width=255,precision=4))
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

    '''
    # データがないときときは解析しない
    if data.size == 0:
        return None

    # サンプル数が１個以下のときは解析しない
    if len(data[0]) <= 1:
        return None

    # 欠損値対策
    data = data[:, ~np.isnan(data).any(0)]

    # 二群の情報を計算
    g_set = np.unique(data[0])
    g1 = data[0] == g_set[0]
    g2 = ~g1
    n1 = np.sum(g1)
    n2 = np.sum(g2)

    if len(g_set) == 2 and min(n1, n2) > 1:
        # do statistical tests
        # note: km should be executed before wilcox and logrank because these
        # tests destruct the order of data
        if km:
            a_stat_km_g1 = stat.km(data[2][g1], data[1][g1])
            a_stat_km_g2 = stat.km(data[2][g2], data[1][g2])
        else:
            a_stat_km_g1, a_stat_km_g2 = None, None
        a_stat_wilcox = stat.gen_wilcox_numpy(data)
        a_stat_logrank = stat.logrank_numpy(data)
        return (a_stat_wilcox, a_stat_logrank, a_stat_km_g1, a_stat_km_g2)
    else:
        return None  # (None,None,None,None)


def gen_survivalTest_with_thretholds(value, event, time, thretholds):
    ''' 閾値を刻みながらvalue vectorを二群に分けて生存分析を行うジェネレータ関数.
    input:
     value,event,time: numpy.array
    output: yield (threshold, n1, n2, results)
    example:
    >>> try:
    ...     from future_builtins import map, zip
    ... except ImportError:
    ...     pass
    ...
    >>> import numpy as np
    >>> from genecoder.lab import analyze
    >>> e = np.array([1,0,1,1,0,0,1,0,0,0,1,1,0,0,0,1,0,0,1,1,0,1,0,0,1,0,1,0,0,0])
    >>> t = np.array([2,20,5,1,3,17,2,3,15,14,12,13,11,11,10,8,8,3,7,3,6,2,5,4,2,3,1,3,2,1])
    >>> v = np.array([float(x)/len(e)/2 for x in range(len(e))])
    >>> for threshold,n1,n2,a_result in analyze.gen_survivalTest_with_thretholds(v,e,t,list(map(
    ... lambda x: float(x)/100, range(5,100,5)))):
    ...     print('[{0}]threshold={1},n1={2},n2={3}'.format('OK' if a_result else 'NOT CALCULATED',
    ... threshold,n1,n2))
    ...
    [OK]threshold=0.05,n1=27,n2=3
    [OK]threshold=0.1,n1=24,n2=6
    [OK]threshold=0.15,n1=21,n2=9
    [OK]threshold=0.2,n1=18,n2=12
    [OK]threshold=0.25,n1=15,n2=15
    [OK]threshold=0.3,n1=12,n2=18
    [OK]threshold=0.35,n1=9,n2=21
    [OK]threshold=0.4,n1=6,n2=24
    [OK]threshold=0.45,n1=3,n2=27
    [NOT CALCULATED]threshold=0.5,n1=0,n2=30
    [NOT CALCULATED]threshold=0.55,n1=0,n2=30
    [NOT CALCULATED]threshold=0.6,n1=0,n2=30
    [NOT CALCULATED]threshold=0.65,n1=0,n2=30
    [NOT CALCULATED]threshold=0.7,n1=0,n2=30
    [NOT CALCULATED]threshold=0.75,n1=0,n2=30
    [NOT CALCULATED]threshold=0.8,n1=0,n2=30
    [NOT CALCULATED]threshold=0.85,n1=0,n2=30
    [NOT CALCULATED]threshold=0.9,n1=0,n2=30
    [NOT CALCULATED]threshold=0.95,n1=0,n2=30

    '''
    for a_threshold in thretholds:
        group = value >= a_threshold
        n1 = np.sum(group)
        n2 = np.sum(~group)
        data = np.array((group, event, time), dtype=float)

        a_result = survivalTest(data)
        yield (a_threshold, n1, n2, a_result)


def analyze_survivalTest(out_file, values, events, times, seq_category, pre_data=OrderedDict(),
                         names=None, thresholds=list(map(lambda x: x / 100, range(5, 100, 5))),
                         drawGraph=False, treatWildTypeAs3rdAxis=False, outDirKM=None):
    ''' テーブルデータ(値が文字列で格納された行列)を数値行列に変換し生存分析を行う関数.
    input: table(values,events,times), seq_category=[mutated/wildtype]
    output: out_file (csv format)

    treatWildTypeAs3rdAxis: WILD TYPEを第3の軸として扱うか. Trueのとき, WILD TYPE以外のデータで検定を行い、
                            WILD TYPEは第3の線として描かれる
    '''
    def gen_event():
        for a_event in events:
            if a_event == 'STOP':
                yield 0
            elif a_event == 'EVENT':
                yield 1
            else:
                yield np.nan

    def gen_array(arr):
        for v in arr:
            if v != '':
                yield float(v)
            else:
                yield np.nan

    # value, time are directly translated from str to float by numpy
    event_np = np.array(list(gen_event()), dtype=np.float64)
    value_np = np.array(list(gen_array(values)), dtype=np.float64)
    time_np = np.array(list(gen_array(times)), dtype=np.float64)
    name_np = np.array(list(range(len(event_np))), dtype=np.float64)

    # サンプル数が１個以下の場合は生存分析を行わない
    if len(event_np) <= 1:
        return False

    # ひとつの大きな２次元表を作成．横方向にサンプルが並ぶ．
    # data[0]: value
    # data[1]: event
    # data[2]: time
    # data[3]: column index
    data = np.array((value_np, event_np, time_np, name_np), dtype=np.float64)

    # 欠損値対策．無効な値が一つでもある列(サンプル)は省く
    data = data[:, ~np.isnan(data).any(0)]

    # WILD TYPE maskを作成
    wildtype_mask = np.array([seq_category[int(index)] == 'wildtype' for index in data[3]])

    # 普通の生存曲線(NonWildType VS WildType)を描画. 各群は２個以上の要素が必須．満たないときは描画しない．
    n1 = np.sum(wildtype_mask)
    n2 = np.sum(~wildtype_mask)
    if outDirKM and n1 > 1 and n2 > 1:
        # do survivalTest
        # input: (group,event,time)
        # output: (stat_wilcox,stat_logrank,stat_g1,stat_g2)
        # TODO: dtype=floatはnp.float64ではだめだった理由を検証
        a_data_standard = np.array((wildtype_mask, data[1], data[2]), dtype=float)
        a_result_standard = survivalTest(a_data_standard)
        # output result
        if a_result_standard:
            (a_stat_wilcox, a_stat_logrank, a_stat_km_g1, a_stat_km_g2) = a_result_standard

            # wildtype
            res_axis1 = stat.km(
                a_data_standard[2, wildtype_mask], a_data_standard[1, wildtype_mask])
            # mutated
            res_axis2 = stat.km(
                a_data_standard[2, ~wildtype_mask], a_data_standard[1, ~wildtype_mask])

            # legend
            legend = ['WildType (n={0})'.format(
                np.sum(wildtype_mask)), 'Mutation (n={0})'.format(np.sum(~wildtype_mask))]
            # pvalue
            pvalue_label = ['$P$ (wilcoxon)={0:1.3}'.format(a_stat_wilcox[3]),
                            '$P$ (logrank)={0:1.3}'.format(a_stat_logrank[1])]
            label = ('months', 'cumulative probability')

            # 5年生存率
            rate = [0.0, 0.0]
            for a_res_axis, a_rate_idx in [(res_axis1, 0), (res_axis2, 1)]:
                pre_t = -1
                for t, a_rate in zip(a_res_axis[0], a_res_axis[3]):
                    if pre_t <= 60 and 60 <= t:
                        rate[a_rate_idx] = a_rate
                        break
                    pre_t = t

            # 中央値
            a_median = [np.median(a_data_standard[2, wildtype_mask]), np.median(
                a_data_standard[2, ~wildtype_mask])]

            a_graph_dir = outDirKM
            a_graph_absfile = os.path.join(
                a_graph_dir,
                'Kaplan-Meier__5year_axis1={0:.3}_axis2={1:.3}__median_axis1={2}_axis2={3}.eps'
                .format(rate[0], rate[1], a_median[0], a_median[1]))
            if not os.path.isdir(a_graph_dir):
                os.makedirs(a_graph_dir)

            stat.plot_km([res_axis1, res_axis2], legend=legend, color=True, label=label, title=(
                '', 10), pvalues=pvalue_label, outfile=a_graph_absfile)

    if treatWildTypeAs3rdAxis:
        # WILD TYPEは別に管理
        data_wildtype = data[:, wildtype_mask]
        # WILD TYPE以外を抽出
        data = data[:, ~wildtype_mask]
    else:
        dummy_wildtype_mask = np.array([False, ] * len(data[0]))
        data_wildtype = data[:, dummy_wildtype_mask]

    header = list(pre_data.keys()) + ['threshold', 'n1', 'n2', 'nwild', 'wildtype_belonging_group',
                                      'median(t)1', 'median(t)2',
                                      'W', 'var(W)', 'z-value', 'wilcox p-value', 'chi-squared',
                                      'logrank p-value', 'group1', 'group2',
                                      'group1 RC', 'group2 RC']
    with open(out_file, 'w') as fout:
        fout.write(','.join(header) + '\n')
        for a_threshold, n1, n2, a_result in gen_survivalTest_with_thretholds(data[0], data[1],
                                                                              data[2], thresholds):
            if verboseWarning and a_result:
                print('[{0}]statistical test with threshold={1} is done.'.format(
                    'OUTPUT' if a_result else 'SKIPPED', a_threshold))
            if a_result:
                '''
                出力データにおけるgroup1, group2, wildtypeの扱いについて
                [1] wildtypeを第3のグラフとしない場合 (treatWildTypeAs3rdAxis == False)
                    group1,group2: RC>=閾値の群, RC< 閾値の群 (どちらかの群がwildtypeをすべて含んでいる)
                    axis1,axis2: wildtypeを含むほうの群, wildtypeを含まないほうの群 wildtypeがないときは
                    group1,group2の順に対応
                    n1,n2: group1の数, group2の数
                    nwild: wildtypeの数
                    belong_group: 'group1' or 'group2' (wildtypeの属する群の名前)
                    注意事項: group1==axis1とは限らないことに注意
                    注意事項: wildtypeな群がどちらのグループに属しているかはmodeで判断できる
                [2] wildtypeを第3のグラフとする場合 (treatWildTypeAs3rdAxis == True)
                    group1,group2: RC>=閾値の群, RC< 閾値の群 (どちらの群もwildtypeを含まない)
                    axis1,axis2,axis3: group1, group2, wildtype
                    n1,n2: group1の数, group2の数
                    nwild: wildtypeの数
                    belong_group: '3rdGroup'
                    注意事項: group1==axis1
                '''
                (a_stat_wilcox, a_stat_logrank,
                 a_stat_km_g1, a_stat_km_g2) = a_result
                g1_mask = (data[0] >= a_threshold)
                g2_mask = ~(data[0] >= a_threshold)
                # nwild = np.sum(wildtype_mask)
                if treatWildTypeAs3rdAxis:
                    belong_group = '3rdGroup'
                else:
                    if wildtype_mask[g1_mask].any():
                        belong_group = 'group1'
                    elif wildtype_mask[g2_mask].any():
                        belong_group = 'group2'
                    else:
                        belong_group = 'noWildType'
                # if not treatWildTypeAs3rdAxis and wildtype_mask[g2_mask].any():
                #    a_temp = (n1,g1_mask,a_stat_km_g1)
                #    (n1,g1_mask,a_stat_km_g1) = (n2,g2_mask,a_stat_km_g2)
                #    (n2,g2_mask,a_stat_km_g2) = a_temp
                # generate output data
                a_line = list(pre_data.values()) + [str(a_threshold), str(n1), str(n2), str(
                    np.sum(wildtype_mask)), belong_group]  # threshold,n1,n2,nwild,belong_group
                a_line += [str(np.median(data[2, g1_mask])),
                           str(np.median(data[2, g2_mask]))]
                a_line += list(map(str, a_stat_wilcox))
                a_line += list(map(str, a_stat_logrank[0:2]))
                a_line += [':'.join([names[int(index)]
                                     for index in data[3, g1_mask]])] if names else ['']
                a_line += [':'.join([names[int(index)]
                                     for index in data[3, g2_mask]])] if names else ['']
                a_line += [':'.join(['{0:.3}'.format(value)
                                     for value in data[0, g1_mask]])] if names else ['']
                a_line += [':'.join(['{0:.3}'.format(value)
                                     for value in data[0, g2_mask]])] if names else ['']
                fout.write(','.join(a_line) + '\n')
                # draw a graph
                if drawGraph:
                    # axis
                    if treatWildTypeAs3rdAxis:
                        res_axis1 = stat.km(data[2, g1_mask], data[1, g1_mask])
                        res_axis2 = stat.km(data[2, g2_mask], data[1, g2_mask])
                        # wild typeがひとつもないときは第3のグラフは描かない
                        if wildtype_mask.any():
                            res_axis3 = stat.km(data_wildtype[2], data_wildtype[1])
                        else:
                            res_axis3 = None
                        legend = ['RC{0}{1} (n={2})'.format(a_axis[0], a_threshold,
                                                            a_axis[1]) for a_axis in (('>=', n1),
                                                                                      ('<', n2))]
                        legend += ['WILD TYPE']
                    else:
                        if belong_group == 'group2':
                            # 変異タイプはTag=Trueとする
                            res_axis2 = stat.km(data[2, g1_mask], data[1, g1_mask],
                                                tag_list=~wildtype_mask[g1_mask])
                            res_axis1 = stat.km(data[2, g2_mask], data[1, g2_mask],
                                                tag_list=~wildtype_mask[g2_mask])
                            legend = ['RC{0}{1} (n={2})'.format(
                                a_axis[0], a_threshold, a_axis[1]) for a_axis in (('<', n2),
                                                                                  ('>=', n1))]
                        else:
                            res_axis1 = stat.km(data[2, g1_mask], data[1, g1_mask],
                                                tag_list=~wildtype_mask[g1_mask])
                            res_axis2 = stat.km(data[2, g2_mask], data[1, g2_mask],
                                                tag_list=~wildtype_mask[g2_mask])
                            legend = ['RC{0}{1} (n={2})'.format(
                                a_axis[0], a_threshold, a_axis[1]) for a_axis in (('>=', n1),
                                                                                  ('<', n2))]
                        res_axis3 = None
                    # 5年生存率
                    rate = [0.0, 0.0]
                    for a_res_axis, a_rate_idx in [(res_axis1, 0), (res_axis2, 1)]:
                        pre_t = -1
                        for t, a_rate in zip(a_res_axis[0], a_res_axis[3]):
                            if pre_t <= 60 and 60 <= t:
                                rate[a_rate_idx] = a_rate
                                break
                            pre_t = t

                    # out dir
                    if treatWildTypeAs3rdAxis:
                        a_graph_dir = os.path.join(
                            os.path.dirname(out_file), 'EXCLUDE_WILDTYPE')
                    else:
                        a_graph_dir = os.path.join(
                            os.path.dirname(out_file), 'INCLUDE_WILDTYPE')
                    s = '{0}_5year_axis1={1:.3}_5year_axis2={2:.3}.eps'.format(a_threshold,
                                                                               rate[0], rate[1])
                    a_graph_absfile = os.path.join(a_graph_dir, s)
                    if not os.path.isdir(a_graph_dir):
                        os.makedirs(a_graph_dir)
                    # pvalue
                    pvalue_label = ['$P$ (wilcoxon)={0:1.3}'.format(a_stat_wilcox[3]),
                                    '$P$ (logrank)={0:1.3}'.format(a_stat_logrank[1])]
                    label = ('months', 'cumulative probability')
                    stat.plot_km([res_axis1, res_axis2, res_axis3], legend=legend, color=True,
                                 label=label, title=('', 10), pvalues=pvalue_label,
                                 outfile=a_graph_absfile)
                    # output the data used in graph
                    for (a_filename, a_data) in [('axis1', data[:, g1_mask]),
                                                 ('axis2', data[:, g2_mask]),
                                                 ('axis3' if res_axis3 else None, data_wildtype)]:
                        if not a_filename:
                            continue
                        a_graph_datafile = os.path.join(
                            a_graph_dir, '{0}_{1}.csv'.format(a_threshold, a_filename))
                        np.savetxt(
                            a_graph_datafile, a_data, delimiter=',', fmt='%f')
                        a_graph_datafile = os.path.join(
                            a_graph_dir, '{0}_{1}_summary.txt'.format(a_threshold, a_filename))

            else:
                # 解析ができなかった場合(二群に分けられなかった場合)は結果を出力しない
                pass

    return True


def analyze_survivalTest_for_database(out_dir, database, target_gf4, target_coders, drawGraph,
                                      databaseStyle='tp53'):
    ''' データベースに対し生存分析を行う関数.
    out_dir: Fasta,RC,Statの各ファイルを格納するディレクトリ
    database: データベースファイル名(カラム群['mutation_id','seq_category','region_name','seq_na',
                                 'RFS(months)','RFS(event)','OS(months)','OS(event)']を含むCSV)
                                 event in {'EVENT','STOP'}, seq_category in {'mutated','wildtype'}
    生存分析用データがあるサンプルのみ計算する
    outFasta: [out_dir]/[basename(database)]/[region]/INPUT.fasta (multiple fasta)
    outRC:    [out_dir]/[basename(database)]/[region]/[GF4]/RC.csv
    outStat:  [out_dir]/[basename(database)]/[region]/[GF4]/[coder]/[event]/STAT.csv
    '''

    # 出力ファイルSTAT.csvのヘッダー
    csv_RC_header = (
        'name',
        'dna sequence',
        'enocoded dna sequence',
        'protein',
        'encoded protein',
        'RC',
        'similarity between original dna sequence and encoded dna sequence',
        'coordinate of GF(4) elements',
        'coder_id',
        'coder_detail')

    # 出力フォルダを微調整
    input_basename = os.path.basename(database)
    out_dir = os.path.join(out_dir, input_basename)
    # 出力フォルダを作成. 既に存在していればそれを削除
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir, ignore_errors=True)
    os.makedirs(out_dir)

    if databaseStyle.lower() == 'tp53':
        # DEFINITIONS
        lifeData_columns = ['mutation_id', 'seq_category', 'region_name',
                            'seq_na', 'RFS(months)', 'RFS(event)', 'OS(months)', 'OS(event)']
        index_mutation_id = 0
        index_seq_category = 1
        index_region_name = 2
        index_seq_na = 3
        index_RFS_months = 4
        index_RFS_event = 5
        index_OS_months = 6
        index_OS_event = 7

        '''
        STEP 1: 生存分析用データを読み込む(高速検索用に辞書としても格納する.
        メモリ効率を考え辞書(key=mutation_id)には同じ行への参照を格納する)
        input: lifeData_columns of database
        output: lifeData, lifeData_dict, region_set, event_set
        {'non_recurrence_months': '', 'event1': '', 'overall_survival_months': '','event2': ''}
        '''
        lifeData_omit = {}
        lifeData = []
        lifeData_dict = OrderedDict()
        with open(database, 'rU') as fin:
            csv_reader = csv.reader(fin)
            for a_line in gen_table(csv_reader=csv_reader, columns=lifeData_columns,
                                    accept_filters=None, omit_filters=lifeData_omit):
                lifeData.append(a_line)
                lifeData_dict[a_line[index_mutation_id]] = a_line
        # 領域の種類を取得
        region_set = set(
            [x[0] for x in gen_columns_from_table(lifeData, [index_region_name])])
        # イベントの種類をセット
        event_set = [('RFS', index_RFS_months, index_RFS_event),
                     ('OS', index_OS_months, index_OS_event)]

        # debug
        print('=== DATABASE SUMMARY ===')
        print('REGIONS: {0}'.format(region_set))
        print('num of samples: {0}'.format(len(lifeData_dict.keys())))
        print('========================')

        # STEP 2: LOOP w.r.t. region
        for a_region in region_set:
            KMgenerated = False
            # 対象の領域の配列のみを選別
            seqs_wrt_region = list(
                gen_columns_from_table(
                    lifeData,
                    [index_mutation_id, index_seq_na],
                    index_value__accept_filters={index_region_name: a_region}
                )
            )

            print('TARGET REGION: {0}, num={1}'.format(
                a_region, len(seqs_wrt_region)))
            if len(seqs_wrt_region) <= 1:
                continue

            # Fastaファイルを出力(確認用). 入力配列はGF4にはよらない共通データであることに注意
            # 出力先: [out_dir]/[basename(database)]_[date]/[region]/INPUT.fasta
            outFasta = os.path.join(out_dir, a_region, 'INPUT.fasta')
            os.mkdir(os.path.dirname(outFasta))
            with open(outFasta, 'w') as fout:
                for a_name, a_seq in seqs_wrt_region:
                    fout.write('>{0}\n{1}\n'.format(a_name, a_seq))
            # STEP 3.1: LOOP w.r.t. GF4 coordinate (このループまででRC距離一覧が計算できる)
            for a_coordinate_of_GF4 in target_gf4:
                # 指定された有限体対応におけるRC距離一覧を計算し出力(その領域の配列数x符号数分).
                # RC距離はイベントによらない共通データであることに注意
                # analyze_RC_distance()と機能的にかぶるが効率重視のためスクラッチで書く.
                # RC距離ファイルのフォーマットを変えるときは注意.
                # 出力先: [out_dir]/[basename(database)]_[date]/[region]/[GF4]/RC.csv
                outRC = os.path.join(
                    out_dir, a_region, a_coordinate_of_GF4, 'RC.csv')
                os.mkdir(os.path.dirname(outRC))
                with open(outRC, 'w') as out_file:
                    csv_writer = csv.writer(out_file)
                    csv_writer.writerow(csv_RC_header)
                    # STEP 3.2: LOOP w.r.t. coder
                    for (coder_id, coder_detail, n, k, a_coder) in map(
                            (lambda x: [x, ] + list(CODERS[x])), target_coders):
                        # 対象の符号によるRC距離を計算(各入力配列における値を逐次計算する)
                        name_RC_dict = OrderedDict()
                        for (a_name, s1, s2, AA1, AA2, RC, simirarity) in gen_RC_distance(
                                seqs=seqs_wrt_region, coder=a_coder,
                                GF4_coordinate=a_coordinate_of_GF4):
                            # RC距離をメモリに格納(必要なメモリは符号数分)
                            name_RC_dict[a_name] = RC
                            # RC距離をファイルに出力(確認用)
                            a_line = [a_name, s1, s2, AA1, AA2, RC, simirarity,
                                      a_coordinate_of_GF4, coder_id, coder_detail]
                            csv_writer.writerow(a_line)
                        # STEP 4: LOOP w.r.t. event
                        # 解析に必要なデータを整理.
                        # name_RC_dict: {配列名:RC距離値}, note that
                        # 配列名は(領域,GF4対応)で絞り込まれている
                        for event_category, event_time_index, event_event_index in event_set:
                            # 統計解析用データを選別
                            # すでに絞り込み済みのname_RC_dict.keys()を利用する
                            stat_value = list(name_RC_dict.values())
                            stat_time = [lifeData_dict[a_name][event_time_index]
                                         for a_name in name_RC_dict.keys()]
                            stat_event = [lifeData_dict[a_name][event_event_index]
                                          for a_name in name_RC_dict.keys()]
                            seq_category = [lifeData_dict[a_name][index_seq_category]
                                            for a_name in name_RC_dict.keys()]
                            for (a_statOutName, a_treatWildType) in [
                                    ('STAT.csv', False),
                                    ('STAT_WILDTYPE_AS_3RD_DATA.csv', True)]:
                                # 統計解析
                                outStat = os.path.join(
                                    out_dir, a_region, a_coordinate_of_GF4, coder_id,
                                    event_category, a_statOutName)
                                # KMグラフ出力先:
                                # [out_dir]/[basename(database)]_[date]/[region]/[event]
                                if not KMgenerated:
                                    outDirKM = os.path.join(
                                        out_dir, a_region, 'Kaplan-Meier WT vs Non WT',
                                        event_category)
                                else:
                                    outDirKM = None
                                if not os.path.isdir(os.path.dirname(outStat)):
                                    os.makedirs(os.path.dirname(outStat))
                                analyze_survivalTest(out_file=outStat, values=stat_value,
                                                     events=stat_event, times=stat_time,
                                                     seq_category=seq_category,
                                                     pre_data=OrderedDict(
                                                         (('gf4', a_coordinate_of_GF4),
                                                          ('coder_id', coder_id))),
                                                     names=list(name_RC_dict.keys()),
                                                     drawGraph=drawGraph,
                                                     treatWildTypeAs3rdAxis=a_treatWildType,
                                                     outDirKM=outDirKM)
                        KMgenerated = True
        # mearge csv
        # 現在, WILDTYPEを取り除いて行った検定に関しては結果を収集せずグラフのみ描画する
        # ここで回収してしまうと, 区別がつかない行ができてしまうため
        # WILDTYPEを取り除いて行った検定の結果に関しては, 最初から取り除いたデータを用意して
        # プログラムを実行すれば得られる
        outStatAll = os.path.join(out_dir, 'STAT_ALL.csv')
        with open(outStatAll, 'w') as out_file:
            header_all = ['region', 'gf4', 'coder', 'event'] + \
                         ['threshold', 'n1', 'n2', 'nwild', 'belong_group',
                          'median(t)1', 'median(t)2', 'W', 'var(W)', 'z-value',
                          'wilcox p-value', 'chi-squared', 'logrank p-value', 'group1',
                          'group2', 'group1 RC', 'group2 RC']
            csv_writer = csv.writer(out_file)
            csv_writer.writerow(header_all)
            for path, dirs, files in sorted(list(os.walk(out_dir))):
                for a_file in files:
                    if a_file != 'STAT.csv':
                        continue
                    # out_dir以降をサブフォルダ名,csvファイル名と分解
                    full_name = os.path.join(path, a_file)
                    name_nodes = full_name[
                        len(os.path.abspath(out_dir)) + len(os.path.sep):].split(os.path.sep)
                    with open(full_name, 'rU') as fin:
                        csv_reader = csv.reader(fin)
                        # ヘッダを空読み
                        next(csv_reader)
                        for a_line in csv_reader:
                            # a_line[gf4,coder_id]はカット
                            csv_writer.writerow(
                                name_nodes[:-1] + a_line[2:])
                    # remove
                    os.remove(full_name)

        # compress and erase small files (graph files, small STAT.csv and
        # directories)
        for a_region in region_set:
            # pack image data w.r.t. a_region
            pack(os.path.join(out_dir, a_region), erase=True)


def fracToStr(frac):
    ''' fraction型を文字列型に変換する関数. '''
    if not FracStyle:
        return str(frac)
    else:
        return str(float(frac))


def gen_RC_distance(seqs, coder, GF4_coordinate='ATGC'):
    ''' 指定した配列群におけるそれぞれの配列のRC距離を計算し返すジェネレータ関数.
    input:
    output: (name, original sequence(length fixed), encoded sequence, AA1, AA2, RC, similarity)
    example:
    >>> import genecoder.lab.analyze
    >>> from genecoder.lab.codec import Coder_Cyclic
    >>> seqs = [('name','atgcatgcatgc')]
    >>> coder = Coder_Cyclic(g_x=[1,1],n=3)
    >>> for a_result in genecoder.lab.analyze.gen_RC_distance(seqs,coder,'ATGC'):
    ... 	print(a_result)
    ...
    ('name', 'atgcatgcatgc', 'ATTCACGCTTGC', 'MHAC', 'IHAC', '0.75', '0.0')

    '''
    bio.set_elements(GF4_coordinate)
    for name1, s1_original in seqs:
        # STEP1: make an information block
        if verboseWarning and len(s1_original) % coder.n != 0:
            print(('[WARNING]omit the tail of {0} '
                   'because the length does not fit with the coder.').format(name1))
        s1 = s1_original[:(len(s1_original) // coder.n) * coder.n]
        if len(s1) == 0:
            if not quietErr:
                print(
                    '[ERROR]the length of {0} is too small for the coder.'.format(name1))
            continue
        s1_codon12 = bio.get_codon12(s1)
        s1_codon12_GF = bio.NA2GF(s1_codon12, table=bio.NA2GF_table)
        # STEP2: encode
        s2_GF = coder.encode(s1_codon12_GF)
        s2 = bio.codon_sort(
            bio.GF2NA(s2_GF, table=bio.GF2NA_table), n=coder.n, k=coder.k)
        # output
        yield (name1, s1, s2, bio.NA2AA(s1), bio.NA2AA(s2), fracToStr(bio.RC(s1, s2)),
               fracToStr(bio.get_similarity(s1, s2)))
