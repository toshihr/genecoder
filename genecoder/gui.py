import os
import itertools
from io import StringIO
from PySide import QtCore, QtGui
from genecoder.mainwindow import Ui_MainWindow
from genecoder.main import coder_list, parameters, analyze_and_output
from genecoder.lab.fasta import Fasta


class MyMainWindow(QtGui.QMainWindow):

    '''
    TODO: Drag&DropでCodeの順番を変えられるようにする
    '''

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
        for a_code in coder_list.keys():
            item = QtGui.QStandardItem(a_code)
            item.setEditable(False)
            item.setToolTip(coder_list[a_code][0])
            if a_code in parameters['coders']:
                self.model_enable.appendRow(item)
            else:
                self.model_disable.appendRow(item)

        # initialize coordinate of GF4 elements
        if len(parameters['GF4']) == 24:
            self.ui.checkBox_GF4_all.setChecked(True)
        else:
            parameters['GF4'] = [parameters['GF4'][-1], ]
            self.ui.lineEdit_GF4.setText(parameters['GF4'][0])

        # initialize filename
        self.ui.lineEdit_output.setText(parameters['output'])

        # initialize input sequences
        self.refreshSequences()

    @QtCore.Slot()
    def run(self):
        # treat coordinate of GF4 elements
        if self.ui.checkBox_GF4_all.isChecked():
            parameters['GF4'] = [''.join(x)
                                 for x in itertools.permutations('ATGC')]
        else:
            parameters['GF4'] = [self.ui.lineEdit_GF4.text().upper(), ]

        # mode 0
        if self.ui.stackedWidget_mode.currentIndex() == 0:
            parameters['mode'] = 'stat'
            if len(parameters['input database']) == 0:
                QtGui.QMessageBox.information(
                    self, 'Error', 'Need to select input database.')
            elif len(parameters['output stat']) == 0:
                QtGui.QMessageBox.information(
                    self, 'Error', 'Need to select output directory.')
            else:
                parameters['output stat'] = os.path.abspath(
                    parameters['output stat'])
                self.ui.lineEdit_output.setText(parameters['output stat'])
                analyze_and_output()
                QtGui.QMessageBox.information(
                    self, 'Finished', 'The results are stored in {0}'
                    .format(parameters['output stat']))

        # mode 1
        else:
            parameters['mode'] = 'standard'
            self.reconstructInputSequences()
            if len(parameters['input sequences']) == 0:
                QtGui.QMessageBox.information(
                    self, 'Error', 'No sequences to analyze.')
            elif len(parameters['output']) == 0:
                QtGui.QMessageBox.information(
                    self, 'Error', 'Need to put output file name.')
            else:
                parameters['output'] = os.path.abspath(parameters['output'])
                self.ui.lineEdit_output.setText(parameters['output'])
                analyze_and_output()
                QtGui.QMessageBox.information(
                    self, 'Finished', 'The results are stored in {0}'.format(parameters['output']))
                # self.ok = True
                # self.close()

    # Fasta Frame
    @QtCore.Slot()
    def clear_sequences(self):
        parameters['input sequences'] = []
        self.refreshSequences()

    @QtCore.Slot()
    def reconstructInputSequences(self):
        fasta = Fasta(StringIO(self.ui.plainTextEdit_input.toPlainText()))
        parameters['input sequences'] = list(fasta.blocks)

    @QtCore.Slot()
    def refreshSequences(self):
        self.ui.plainTextEdit_input.clear()
        for (a_name, a_seq) in parameters['input sequences']:
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
            parameters['coders'].append(a_item.text())

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
            parameters['coders'].remove(a_item.text())

    # mode 0
    @QtCore.Slot()
    def selectInputDatabase(self):
        file_name = QtGui.QFileDialog.getOpenFileName(
            self, self.tr('Select Input File'), '', self.tr('CSV (*.csv);;any file (*.*)'))[0]
        if file_name != '':
            parameters['input database'] = file_name
            self.ui.lineEdit_inputDatabase.setText(
                parameters['input database'])

    @QtCore.Slot()
    def selectOutDir(self):
        self.reconstructInputSequences()
        dir_name = QtGui.QFileDialog.getExistingDirectory(
            self, self.tr('Select Output Directory'))
        if dir_name != '':
            parameters['output stat'] = dir_name
            self.ui.lineEdit_output_stat.setText(parameters['output stat'])

    @QtCore.Slot()
    def outDirChanged(self, text):
        parameters['output stat'] = text

    @QtCore.Slot()
    def inputDatabaseChanged(self, text):
        parameters['input database'] = text

    # mode 1
    @QtCore.Slot()
    def selectFastaFile(self):
        self.reconstructInputSequences()
        file_name = QtGui.QFileDialog.getOpenFileName(self, self.tr(
            'Select Input File'), '', self.tr('FASTA (*.fasta);;any file (*.*)'))[0]
        if file_name != '':
            fasta = Fasta(open(file_name))
            self.ui.logTextBrowser.append(str(fasta.blocks))
            parameters['input sequences'] += fasta.blocks
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
            parameters['output'] = file_name
            self.ui.lineEdit_output.setText(parameters['output'])

    @QtCore.Slot()
    def outFileChanged(self, text):
        parameters['output'] = text
