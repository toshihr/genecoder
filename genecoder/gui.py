# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals
import os
import itertools
from io import StringIO
from PySide import QtCore, QtGui
from genecoder.mainwindow import Ui_MainWindow
from genecoder.lab.fasta import Fasta
from genecoder.resource import CODERS
from genecoder.main import mode_stat, mode_distance

pipe_args = {}


class MyMainWindow(QtGui.QMainWindow):

    def __init__(self, parent=None):
        super(MyMainWindow, self).__init__(parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # resize
        self.resize(1024, 750)

        # disable log
        self.ui.logDock.hide()

        # connect
        # mode selection
        QtCore.QObject.connect(self.ui.buttonGroup, QtCore.SIGNAL(
            "buttonClicked(int)"), self.ui.stackedWidget_mode.setCurrentIndex)
        # NOTE: mainwindow.pyを自動生成した場合、下記該当個所の修正が必要. IDを1,2とする
        # self.buttonGroup.addButton(self.radioButton_mode0,0)
        # self.buttonGroup.addButton(self.radioButton_mode1,1)
        # coder selection
        QtCore.QObject.connect(
            self.ui.pushButton_toEnable, QtCore.SIGNAL('clicked()'), self.toEnable)
        QtCore.QObject.connect(
            self.ui.pushButton_toDisable, QtCore.SIGNAL('clicked()'), self.toDisable)
        # mode 0 (Survival Analysis)
        QtCore.QObject.connect(self.ui.toolButton_open_inputDatabase, QtCore.SIGNAL(
            'clicked()'), self.selectInputDatabase)
        QtCore.QObject.connect(
            self.ui.toolButton_open_outdir, QtCore.SIGNAL('clicked()'), self.selectOutDir)
        QtCore.QObject.connect(self.ui.lineEdit_inputDatabase, QtCore.SIGNAL(
            'textChanged(QString)'), self.inputDatabaseChanged)
        QtCore.QObject.connect(self.ui.lineEdit_output_stat, QtCore.SIGNAL(
            'textChanged(QString)'), self.outDirChanged)
        # mode 1 (Code Analysis)
        QtCore.QObject.connect(
            self.ui.toolButton_open, QtCore.SIGNAL('clicked()'), self.selectFastaFile)
        QtCore.QObject.connect(
            self.ui.toolButton_clear, QtCore.SIGNAL('clicked()'), self.clear_sequences)
        QtCore.QObject.connect(
            self.ui.toolButton_open_outfile, QtCore.SIGNAL('clicked()'), self.selectOutFile)
        QtCore.QObject.connect(self.ui.lineEdit_output, QtCore.SIGNAL(
            'textChanged(QString)'), self.outFileChanged)
        # run button
        QtCore.QObject.connect(
            self.ui.pushButton_run, QtCore.SIGNAL('clicked()'), self.run)

        # initialize codeListView
        # inner containers
        self.model_disable = QtGui.QStandardItemModel()
        self.model_enable = QtGui.QStandardItemModel()
        # GUI
        self.ui.codeListView_disable.setModel(self.model_disable)
        self.ui.codeListView_enable.setModel(self.model_enable)
        self.ui.codeListView_disable.setSelectionMode(
            QtGui.QAbstractItemView.ContiguousSelection)
        self.ui.codeListView_enable.setSelectionMode(
            QtGui.QAbstractItemView.ContiguousSelection)

        # initialize codeListView from parameters['coders']
        for a_code in CODERS.keys():
            item = QtGui.QStandardItem(a_code)
            item.setEditable(False)
            item.setToolTip(CODERS[a_code][0])
            if a_code in pipe_args['--coder']:
                self.model_enable.appendRow(item)
            else:
                self.model_disable.appendRow(item)

        # initialize coordinate of GF4 elements
        if len(pipe_args['--gf4']) == 24:
            self.ui.checkBox_GF4_all.setChecked(True)
        else:
            if(len(pipe_args['--gf4']) > 0):
                self.ui.lineEdit_GF4.setText(pipe_args['--gf4'][0])

        # initialize filename
        if pipe_args['distance']:
            self.ui.lineEdit_output.setText(pipe_args['--output'])
        elif pipe_args['stat']:
            self.ui.lineEdit_output.setText(pipe_args['--outdir'])

        # initialize input sequences
        self.refreshSequences()

    @QtCore.Slot()
    def run(self):
        # treat coordinate of GF4 elements
        if self.ui.checkBox_GF4_all.isChecked():
            pipe_args['--gf4'] = [''.join(x) for x in itertools.permutations('ATGC')]
        else:
            pipe_args['--gf4'] = [self.ui.lineEdit_GF4.text().upper(), ]

        # mode 0 (stat mode)
        if self.ui.stackedWidget_mode.currentIndex() == 0:
            pipe_args['distance'] = None
            pipe_args['stat'] = True
            if len(pipe_args['--input']) == 0:
                QtGui.QMessageBox.information(self, 'Error', 'Need to select input database.')
            elif len(pipe_args['--outdir']) == 0:
                QtGui.QMessageBox.information(self, 'Error', 'Need to select output directory.')
            else:
                pipe_args['--ourdir'] = os.path.abspath(pipe_args['--ourdir'])
                self.ui.lineEdit_output.setText(pipe_args['--ourdir'])

                mode_stat(pipe_args)

                QtGui.QMessageBox.information(self, 'Finished', 'The results are stored in {0}'
                                                    .format(pipe_args['--outdir']))

        # mode 1 (distance mode)
        else:
            pipe_args['distance'] = True
            pipe_args['stat'] = None
            self.reconstructInputSequences()
            if len(pipe_args['--input']) == 0:
                QtGui.QMessageBox.information(self, 'Error', 'No sequences to analyze.')
            elif len(pipe_args['--output']) == 0:
                QtGui.QMessageBox.information(self, 'Error', 'Need to put output file name.')
            else:
                pipe_args['--output'] = os.path.abspath(pipe_args['--output'])
                self.ui.lineEdit_output.setText(pipe_args['--output'])

                mode_distance(pipe_args)

                QtGui.QMessageBox.information(
                    self, 'Finished', 'The results are stored in {0}'.format(pipe_args['--output']))
                # self.ok = True
                # self.close()

    # Fasta Frame
    @QtCore.Slot()
    def clear_sequences(self):
        pipe_args['input sequences'] = []
        self.refreshSequences()

    @QtCore.Slot()
    def reconstructInputSequences(self):
        fasta = Fasta(StringIO(self.ui.plainTextEdit_input.toPlainText()))
        pipe_args['input sequences'] = list(fasta.blocks)

    @QtCore.Slot()
    def refreshSequences(self):
        self.ui.plainTextEdit_input.clear()
        for (a_name, a_seq) in pipe_args['input sequences']:
            self.ui.plainTextEdit_input.appendPlainText('>' + a_name)
            self.ui.plainTextEdit_input.appendPlainText(a_seq)

    @QtCore.Slot()
    def toEnable(self):
        # get the selected item indexes
        indexes = self.ui.codeListView_disable.selectedIndexes()
        # move the selected items from "self.model_disable" container to
        # "items" list
        items = []
        for i in indexes[::-1]:
            items.insert(0, self.model_disable.takeItem(i.row()))
            self.model_disable.removeRow(i.row())
        # move the selected items from the "items" list to "self.model_enable"
        # container
        for a_item in items:
            self.model_enable.appendRow(a_item)
            pipe_args['--coder'].append(a_item.text())

    @QtCore.Slot()
    def toDisable(self):
        # get the selected item indexes
        indexes = self.ui.codeListView_enable.selectedIndexes()
        # move the selected items from "self.model_enable" container to "items"
        # list
        items = []
        for i in indexes[::-1]:
            items.insert(0, self.model_enable.takeItem(i.row()))
            self.model_enable.removeRow(i.row())
        # move the selected items from the "items" list to "self.model_disable"
        # container
        for a_item in items:
            self.model_disable.appendRow(a_item)
            pipe_args['--coder'].remove(a_item.text())

    # mode 0
    @QtCore.Slot()
    def selectInputDatabase(self):
        file_name = QtGui.QFileDialog.getOpenFileName(
            self, self.tr('Select Input File'), '', self.tr('CSV (*.csv);;any file (*.*)'))[0]
        if file_name != '':
            pipe_args['--input'] = file_name
            self.ui.lineEdit_inputDatabase.setText(pipe_args['--input'])

    @QtCore.Slot()
    def selectOutDir(self):
        self.reconstructInputSequences()
        dir_name = QtGui.QFileDialog.getExistingDirectory(
            self, self.tr('Select Output Directory'))
        if dir_name != '':
            pipe_args['--ourdir'] = dir_name
            self.ui.lineEdit_output_stat.setText(pipe_args['--ourdir'])

    @QtCore.Slot()
    def outDirChanged(self, text):
        pipe_args['--ourdir'] = text

    @QtCore.Slot()
    def inputDatabaseChanged(self, text):
        pipe_args['--input'] = text

    # mode 1
    @QtCore.Slot()
    def selectFastaFile(self):
        self.reconstructInputSequences()
        file_name = QtGui.QFileDialog.getOpenFileName(self, self.tr(
            'Select Input File'), '', self.tr('FASTA (*.fasta);;any file (*.*)'))[0]
        if file_name != '':
            fasta = Fasta(open(file_name))
            self.ui.logTextBrowser.append(str(fasta.blocks))
            pipe_args['input sequences'] += [fasta.blocks]
            self.refreshSequences()
            self.ui.logTextBrowser.append(file_name)

    @QtCore.Slot()
    def selectOutFile(self):
        self.reconstructInputSequences()
        file_name = QtGui.QFileDialog.getSaveFileName(self, self.tr(
            'Select Output File'), '', self.tr('CSV (*.csv);;any file (*.*)'))[0]
        if file_name != '':
            if file_name[-4:] != '.csv':
                file_name += '.csv'
            pipe_args['--output'] = file_name
            self.ui.lineEdit_output.setText(pipe_args['--output'])

    @QtCore.Slot()
    def outFileChanged(self, text):
        pipe_args['--output'] = text
