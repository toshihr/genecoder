==========================================================
genecoder: Code analyzer for the coding region of a gene.
==========================================================

.. image:: https://travis-ci.org/kerug/genecoder.svg
    :target: https://travis-ci.org/kerug/genecoder


Requirements
============

- Qt4

Qt4 (Mac OS X)
--------------

Install Qt4 via Homebrew is recommended. Install Homebrew_ is as follows:

::

    $ ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

.. _Homebrew: http://brew.sh/


Install Qt4 is as follows:

::

    $ brew install qt

Qt4 (Ubuntu)
------------

::

    $ sudo apt-get install libqt4-dev

How to install
==============

The genecoder can be installed via pip_.

.. _pip: https://pip.pypa.io/en/latest/installing.html

::

    $ pip install genecoder
    $ pyside_postinstall.py -install

Usage examples
==============

Calculate RC distance:


::

    $ genecoder distance --coder n3_1 --coder n3_2 -gf4 ATGC --seq label1:ATGCATGCATGC --output [result]
    $ genecoder --all --gf4all --input [FASTA]

Survival analysis:

::

    $ genecoder stat --graph --coder n3_1 --outdir [result_dir] --input [tp53 database file]


The results are stored in ``result_dir`` folder.


Generate FASTA file from csv database:

::

    $ genecoder csv2fasta <idx_name> <idx_seq> [<length>] [--input=<csv>] [--output=<output>]

Use GUI:

::

    $ genecoder gui

Show help:

::

    $ genecoder -h

Show support coders:

::

    $ genecoder list


TP53 database file format
=========================

TP53 database is a CSV(Comma-Separated Values) format file.
Columns should be the followings:

- mutation_id
- seq_category
- region_name
- seq_na
- RFS(months)
- RFS(event)
- OS(months)
- OS(event)


How to develop
==============

Developers should use pyenv_ and `pyenv-virtualenv`_.

.. _pyenv: https://github.com/yyuu/pyenv
.. _pyenv-virtualenv: https://github.com/yyuu/pyenv-virtualenv

Mac OS X users can be installed via Homebrew:

::

    $ brew install pyenv-virtualenv


How to construct an environment:

::

    $ git clone https://github.com/kerug/genecoder.git
    $ cd genecoder

    $ pyenv install 2.7.5
    $ pyenv install 3.4.1
    $ pyenv virtualenv 2.7.5 genecoder-2.7.5
    $ pyenv virtualenv 3.4.1 genecoder-3.4.1
    $ pyenv local genecoder-2.7.5 genecoder-3.4.1

    $ pip install -r test-requirements.txt
    $ pip3 install -r test-requirements.txt


Tests for Python 2 & 3:

::

    $ tox


Alternatively,

::

    $ python setup.py test


Sometimes, the following commands are needed:

::

    $ pyside_postinstall.py -install
    $ pyenv rehash


Qt creator's user-interface (\*.ui) can be converted to python code as follows:

::

    $ pyside-uic -o mainwindow.py mainwindow.ui


References
==========

- Sato Keiko, Toshihide Hara, and Masanori Ohya. "The code structure of the p53 DNA-binding domain
  and the prognosis of breast cancer patients." Bioinformatics 29.22 (2013): 2822-2825. [Link_]
- http://theory.cs.uvic.ca/gen/poly.html

.. _Link: http://www.ncbi.nlm.nih.gov/pubmed/23986567
