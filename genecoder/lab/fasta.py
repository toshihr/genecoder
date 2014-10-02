from __future__ import absolute_import, division, print_function, unicode_literals
import subprocess


def msf2fasta(inFile, outFile):
    args = ['/usr/bin/seqret', '-auto', inFile, 'fasta:' + outFile]
    try:
        p = subprocess.Popen(
            args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except OSError:
        print('Failed to execute command: ' + args[0])
    (stdoutdata, stderrdata) = p.communicate()


class Fasta:

    def __init__(self, inFileIO, aligned=False):
        self.inFileIO = inFileIO
        self.blocks = []
        self.num = 0
        self.id_global = 0.0
        self.id_local = 0.0
        self.read()
        if aligned:
            self.calc_id()

    def read(self):
        self.blocks = []
        self.num = 0
        seq_buffer = []
        name = ''
        for line in self.inFileIO:
            # skip empty line
            if line is None:
                continue
            if len(line) == 0:
                continue
            # read
            if line[0] == '>':
                if name != '':
                    self.blocks.append((name, ''.join(seq_buffer)))
                # start new fasta block
                name = line[1:].rstrip()
                self.num += 1
                seq_buffer = []
            else:
                seq_buffer.append(line.rstrip())
        if name != '':
            self.blocks.append((name, ''.join(seq_buffer)))

    def calc_id(self):
        match = 0
        num_gap = 0
        length = len(self.blocks[0][1])
        for index in range(0, length):
            site_gapped = False
            site_match = True
            symbol = ''
            for name, seq in self.blocks:
                c = seq[index]
                if c == '-':
                    num_gap += 1
                    site_gapped = True
                    break
                if symbol == '':
                    symbol = c
                elif symbol != c:
                    site_match = False
                    break
            if not site_gapped and site_match:
                match += 1
        self.id_global = match / length
        self.id_local = match / (length - num_gap)
