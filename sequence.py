__author__ = "asyavuz, sozerberk"
#
# sequence.py
# Sequence wrapper class for handling sequences easily, and calculating sequence-level features
#
# Authors: Ahmet Sinan Yavuz <asinanyavuz@sabanciuniv.edu>, Namik Berk Sozer <namikberk.sozer@std.yeditepe.edu.tr>
#
# Copyright (c) 2015 Ahmet Sinan Yavuz, Namik Berk Sozer, Osman Ugur Sezerman
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
import re
import residue

# Normal Amino Acid Alphabet
aa_dict = {"ALA": 0, "ARG": 1, "ASN": 2, "ASP": 3, "CYS": 4, "GLN": 5, "GLU": 6, "GLY": 7, "HIS": 8,
           "ILE": 9, "LEU": 10, "LYS": 11, "MET": 12, "PHE": 13, "PRO": 14, "SER": 15, "THR": 16, "TRP": 17,
           "TYR": 18, "VAL": 19}

class Seq(object):
    '''
        Sequence wrapper class.
    '''
    identifier = ''
    sequence = ''
    residues = []
    features = {}

    def __init__(self, identifier, sequence, features=None):
        '''
        Initializes a sequence object
        :param identifier: Sequence identifier
        :param sequence: Amino acid sequence
        :param features: Pre-calculated features
        :return:
        '''
        self.identifier = identifier
        self.sequence = sequence
        self.residues = []
        self.residues.append(None)

        for i in xrange(0, len(self.sequence)): # Generate a residue object for each residue
            tmp = residue.Residue(no=(i + 1), aa=self.sequence[i], parent=self)
            tmp.features["class"] = False
            self.residues.append(tmp)
        if features:
            self.features = features
        else:
            self.features = {}

    def add_iupred(self, filename=''):
        '''
        IUPred results parser
        :param filename: IUPred prediction result file
        :return: Nothing
        '''
        ifh = open(filename, "r")
        lines = ifh.readlines()
        ifh.close()

        iuline = re.compile(r'(\d+)\s*(\w+)\s*(\d+\.\d+)')

        for line in lines:
            ad = iuline.search(line)
            if (ad):
                if self.residues[int(ad.groups(1)[0])].aa == ad.groups(1)[1]:
                    self.residues[int(ad.groups(1)[0])].features['iupred'] = float(ad.groups(1)[2])
                else:
                    raise Exception('Mismatch in IUPRED File with Sequence!')

    def add_pssm(self, filename=''):
        '''
        PSSM  parser
        :param filename: IUPred prediction result file
        :return: Nothing
        '''
        ifh = open(filename, "r")
        lines = ifh.readlines()
        ifh.close()

        aa_order = []

        line_pattern = re.compile(r'\s*(\d+)\s*([A-Z])\s*(.+)?[\n\r]+')
        aa_order_pattern = re.compile(r'\s+([A-Z\s]+)?[\n\r]+')

        for idx, line in enumerate(lines):
            if line == "" or line.startswith("Last position-specific"):
                continue

            aa_order_match = aa_order_pattern.search(line)
            if aa_order_match:
                aa_order = aa_order_match.group(1).split()[0:20]
            else:
                line_match = line_pattern.search(line)
                if line_match:
                    rc = int(line_match.group(1))
                    aa = line_match.group(2)
                    if self.residues[rc].aa != aa:
                        raise Exception("Invalid PSSM File! Amino acid mismatch at line %d for aminoacid %d in file %s.\nHave "
                                      "%s, got %s." % (idx, rc, filename, self.residues[rc].aa, aa))

                    elems = line_match.group(3).split()[0:20]
                    self.residues[rc].features['PSSM'] = [ int(x) for x in elems ]
                else:
                    continue

        self.features['PSSM_AA_Order'] = aa_order

    def set_aa_freq(self):
        '''
        Calculates frequency of each amino acid and amino acid groups in the protein
        :return: Nothing
        '''
        rawamino_acids = ['A', 'N', 'C', 'Q', 'H', 'L', 'M', 'P', 'T', 'Y', 'R', 'D', 'E', 'G', 'I', 'K', 'F', 'S', 'W', 'V']
        rawamino_acids_sz = ['[IVLM]', '[RKH]', '[DE]', '[QN]', '[ST]', '[A]', '[G]', '[W]', '[C]', '[YF]', '[P]']

        attribute_names = ['freqSeq'+i for i in rawamino_acids]
        attribute_vals = [(self.sequence.count(val))/float(len(self.sequence)) for val in rawamino_acids]

        attribute_names.extend(['freqSeq'+i for i in rawamino_acids_sz])
        attribute_vals.extend([len(re.findall(val, self.sequence))/float(len(self.sequence)) for val in rawamino_acids_sz])

        for idx, val in enumerate(attribute_names):
            self.features[val] = attribute_vals[idx]

    def get_fasta(self):
        '''
        Returns FASTA (Pearson) formatted string of the sequence.
        :return: FASTA formatted string
        '''
        fasta = ">%s\n" % (self.identifier)
        seq = self.sequence
        for i in range(0, len(seq), 60):
            fasta += (seq[i:i + 60] + "\n")
        return fasta

    def __getitem__(self, item):
        return self.residues[item]

    def __len__(self):
        return len(self.residues)

    def __iter__(self):
          return iter(self.residues)