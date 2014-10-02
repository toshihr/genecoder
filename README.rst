GeneCoder [WIP]
=========

work in progress.


Requirements
------------

- python & pip
- Qt4


Install Requirements
--------------------

Qt4
~~~

Install Qt4 via Homebrew is recommended. Install Homebrew is as follows:

::

    $ ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

Install Qt4 is as follows:

::

    $ brew install qt

python & pip
~~~~~~~~~~~~


Install
-------

::

    $ pip install git+https://github.com/kerug/genecoder.git


Build
-----

::

    $ pip install -r requirements.txt


File Format
-----------
統計解析を行うための入力ファイルは以下のフォーマット
'mutation_id', 'seq_category', 'region_name',
'seq_na', 'RFS(months)', 'RFS(event)', 'OS(months)', 'OS(event)'
というカラムが存在するCSV(delimiter=',', quotechar='"')ファイル

# ============================================
# - 遺伝子配列(CDS)に対し符号解析を行うソフト -
# MODE default:
#  - coding sequenceにおける情報部(コドン単位における第1,2番目の塩基)を元に符号化を行い、
#    元の配列とのRC距離(翻訳後アミノ酸配列における一致率)を計算する．
#  - 結果はデフォルトでは標準出力に出力される．
#  - 結果はCSV(del=',', quot='"')フォーマット．
#  ex) >./genecoder.py -seq seq1 ATGATGATGATGATG -code bch_n3_1 -code bch_n3_2 -c Iwadare -c SelfOrtho -gf4 ATGC
#  ex) >./genecoder.py -all -gf4all [SEQUENCE FILE (MULTIPLE FASTA)]
#
# MODE -stat, -stat_with_graph:
#  - 変異情報,イベント情報を含むデータを統計解析する．各符号，各RC閾値ごとに二群にわけ，有意差があるものを出力する．
#  - 入力CSVからは下記のカラムが読まれる(位置は関係ない)．これ以外のカラムは無視される．
#    ['mutation_id','region_name','seq_na','RFS(months)','RFS(event)','OS(months)','OS(event)']
#  ex) >./genecoder.py -gf4all -all -stat -quiet -quieterr -output [OUTPUT DIR] [INPUT CSV(del=',', quot='"')]
# MODE -showcoders:
#   現在対応している符号を表示する
# --- main features ---
# 対応する入力ファイル形式: mutiple FASTA
# 対応する符号をプラグインとして追加できる: genecoder.iniにて記述する
# GUIモード
# Linux,Windowsでの動作
# --- compile ---
# >python3 setup.py build_ext --inplace
# --- make an install package ---
# >python3 setup.py build
# >python3 setup.py bdist_msi
# ============================================
#
# === KNOWN BUG ===
# 入力ファイルがないときでもGUIモードが起動できるようにする
#
# === TODO ===
# 2012/9/24現在stat modeがおかしいlife dataのないサンプルに対しても計算に含めてしまう
# stat modeの閾値1つ版をつくる. ついでにgraphを出力するかのオプションをつける.
# graph modeを作る. graph modeのGUIを搭載する. stat modeではdatabase名を選ぶようにしoutputはdirectoryを選ぶようにする. 配列入力はdisableにする.
# パリティビットの挿入法(シャッフル)を選べるようにする
# DNAの長さが符号長とずれている場合に対処するオプションをつける
#
# === MEMO ===
# ::Qt creatorで作成したインターフェイスファイルmainwindow.uiをpython用コードに変更するコマンド(pyside版)
# pyside-uic -o ../mainwindow.py mainwindow.ui
#
# ::原始多項式を求められるサイト
#http://theory.cs.uvic.ca/gen/poly.html
#
# ::数学的内容のメモ
# 符号理論の基礎
# 情報系列 i = (i0,i1,...,ik-1) より対応する多項式はk-1次
# 符号多項式X(t)はn-1次以下の多項式でX(t)=Q(t)G(t)
# 生成多項式G(t)の次数は検査記号数に等しい
#
# ::気になること
# P1. {ATGC} :-> GF(4) の対応にはバリエーションがある
# P2. shuffleも一種の符号化(正確にはソート). これにもバリエーションがある
# (n,k,g(x),P1,P2)の組に対しRCが１つ定まる. 既存研究ではP1,P2は固定している.
# 疑問１：遺伝暗号表的な生体符号はあくまでtandem　データの表現に幅を持たせる符号
# 線形符号はエラー訂正的な符号
# 大事なデータを保護するための戦略が違う。これらを比べることに意味があるのか？
# 疑問２：符号語が近ければ，本来備える符号の仕組みに近いとはいえない。
# 符号語を比較することはナンセンスに感じる。なにか違う比較の指標はないものか。たんちょう度みたいのとか。
# 遺伝暗号表におけるたんちょう度 <-> 人工符号における符号長nに対する情報長kの割合。～速度？　が対応して比較できるのでは
# 今の手法だと符号化速度k/n = 2/3と固定されている。2/3以外の人工符号で行う方法はないか？
# 見方を変えるとアミノ酸にエンコードしたあとのたんちょう度を見るといいのかも。それが遺伝暗号表におけるそれと近いもので比較するとか。
# 疑問３：そもそもCDSをエラー訂正符号と考えるのは無意味では？CDS上にはあくまで冗長性がみられるだけ。
# エラー訂正は5'-3'鎖に対する3'-5'が担っている。どっかのサイトでDNAのエラー訂正に関する記事を読んだ。
# どこかのサイトではイントロン、スプライシングは畳み込み符号的だととらえていた。
# DNA上で実際行われているエラー訂正符号を考え、実際に符号の数理を構築してみる。そのメリットを考える。エラーよりスピード重視とか。
#
# ::検証法の提案
# あくまで(3,2)符号だけにしぼる。こうすることで上記のたんちょう度比較が行える。P1の対応を変えることで
# 検証する符号のバリエーションは増やせる。
# 感じる違和感の多くは(6,4)符号などにしたときのパリティビットの扱い、ソート依存なところにある。それを緩和できる。
