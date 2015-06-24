__author__ = "asyavuz, sozerberk"
#
# residue.py
# Residue wrapper class for calculating window-level features
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
import numpy
import re
import math
import itertools

sz_grouping = { "A": "A", "R": "RKH", "N": "QN", "D": "DE", "C": "C", "E": "DE", "Q": "QN", "G": "G",
            "H": "RKH", "I": "IVLM", "L": "IVLM", "K": "RKH", "M": "IVLM", "F": "YF", "P": "P", "S": "ST",
            "T": "ST", "W": "W", "Y": "YF", "V": "IVLM"}

class Bunch(dict):
    """
    Container object for datasets: dictionary-like object that
    exposes its keys as attributes.
    """

    def __init__(self, **kwargs):
        dict.__init__(self, kwargs)
        self.__dict__ = self


hopp_hydrophobicity = { 'A': -0.50, 'C': -1.00, 'D': 3.00, 'E': 3.00, 'F': -2.50, 'G': 0.00, 'H': -0.50,
                        'I': -1.80, 'K': 3.00, 'L': -1.80, 'M': -1.30, 'N': 0.20, 'P': 0.00, 'Q': 0.20,
                        'R': 3.00,  'S': 0.30, 'T': -0.40, 'V': -1.50, 'W': -3.40, 'Y': -2.30}
class Residue:
    '''
    Residue class for protein sequences. This class is used to calculate window features from central lysine residue.
    '''
    parent = None   # Parent sequence
    no = 0   # Amino acid no in sequence
    aa = ''   # Amino acid
    features = {}   # Unprocessed feature dictionary
    pro_features = None   # Processed feature pointer

    def __init__(self, no=0, aa='', parent=None, raw_features=None):
        '''
        Initiates a residue object
        :param no: Amino acid no
        :param aa: Amino acid
        :param parent: Parent protein sequence
        :param raw_features: Unprocessed feature dictionary
        :return: Nothing
        '''
        self.parent = parent
        self.no = no
        self.aa = aa
        if raw_features:
            self.features = raw_features
        else:
            self.features = {}
        self.pro_features = {}

    def set_selected(self, window_size=21):
        '''
        Wrapper function for calling all feature calculation methods.
        :param window_size: Size of the sequence window
        :return: Dictionary containing feature set short codes as keys and the names of calculated features as values.
        '''
        comp1_attr = self.set_compositions(grouping="no", k=1, window_size=window_size)
        comp1_attr_sz = self.set_compositions(grouping="sezerman", k=1, window_size=window_size)
        iupr_attr = self.set_iupred_real(window_size=window_size)
        iupb_attr = self.set_iupred_binary(window_size=window_size)
        termini_attr = self.set_termini()
        winfo_attr = self.set_hydrophobicity(window_size=window_size)
        aacounts_attr = self.set_aacounts(grouping="no", window_size=window_size)
        aacounts_sz_attr = self.set_aacounts(grouping="sezerman", window_size=window_size)
        aafreqs_attr = self.set_aafreqratios(grouping="no", window_size=window_size)
        aafreqs_sz_attr = self.set_aafreqratios(grouping="sezerman", window_size=window_size)
        pssm_attr = self.set_pssm(window_size=window_size)

        return {'COMP1': comp1_attr,
                'COMP1_SZ': comp1_attr_sz, 'PSSM': pssm_attr,
                'AACNT': aacounts_attr, 'AACNTSZ': aacounts_sz_attr,
                'AAFRQ': aafreqs_attr, 'AAFRQSZ': aafreqs_sz_attr,
                'IUPR': iupr_attr, 'IUPB': iupb_attr,
                'TERMINI': termini_attr,
                'WINFO': winfo_attr
        }

    def set_aacounts(self, grouping="no", window_size=21):
        '''
        Calculates the amount of each amino acid in a window.
        :param grouping: Amino acid groping (no for no grouping, 'sezerman' for sezerman grouping
        :param window_size: Size of the window
        :return: Names of the calculated features
        '''
        rawamino_acids = ['A', 'N', 'C', 'Q', 'H', 'L', 'M', 'P', 'T', 'Y', 'R', 'D', 'E', 'G', 'I', 'K', 'F', 'S', 'W', 'V']
        rawamino_acids_sz = ['[IVLM]', '[RKH]', '[DE]', '[QN]', '[ST]', '[A]', '[G]', '[W]', '[C]', '[YF]', '[P]']

        left, right = self.get_window(window_size)
        window = left+right
        if grouping == "no":
            attribute_names = ['wc'+i for i in rawamino_acids]
            attribute_vals = [window.count(val) for val in rawamino_acids]
        elif grouping == "sezerman":
            attribute_names = (['wc'+i for i in rawamino_acids_sz])
            attribute_vals = [len(re.findall(val, window)) for val in rawamino_acids_sz]

        for idx, val in enumerate(attribute_names):
            self.pro_features[val] = attribute_vals[idx]

        return attribute_names

    def set_aafreqratios(self, grouping="no", window_size=21):
        '''
        Calculates the frequency ratios of each amino acid in a window.
        :param grouping: Amino acid groping (no for no grouping, 'sezerman' for sezerman grouping
        :param window_size: Size of the window
        :return: Names of the calculated features
        '''
        rawamino_acids = ['A', 'N', 'C', 'Q', 'H', 'L', 'M', 'P', 'T', 'Y', 'R', 'D', 'E', 'G', 'I', 'K', 'F', 'S', 'W', 'V']
        rawamino_acids_sz = ['[IVLM]', '[RKH]', '[DE]', '[QN]', '[ST]', '[A]', '[G]', '[W]', '[C]', '[YF]', '[P]']

        left, right = self.get_window(window_size)
        window = left+right

        if grouping == "no":
            attribute_names = ['fr'+i for i in rawamino_acids]
            attribute_vals = []
            for val in rawamino_acids:
                if self.parent.features['freqSeq'+val] == 0:
                    attribute_vals.append(0)
                else:
                    attribute_vals.append((window.count(val)/float(window_size))/self.parent.features['freqSeq'+val])
        elif grouping == "sezerman":
            attribute_names = (['fr'+i for i in rawamino_acids_sz])
            attribute_vals = []
            for val in rawamino_acids_sz:
                if self.parent.features['freqSeq'+val] == 0:
                    attribute_vals.append(0)
                else:
                    attribute_vals.append((len(re.findall(val, window))/float(window_size))/self.parent.features['freqSeq'+val])

        for idx, val in enumerate(attribute_names):
            self.pro_features[val] = attribute_vals[idx]

        return attribute_names

    def set_pssm(self, window_size=7):
        '''
        Adds calculated PSSM values to the window features
        :param window_size: Size of the window
        :return: Names of the calculated features
        '''
        attribute_names = []
        attribute_vals = []
        aa_order = self.parent.features["PSSM_AA_Order"]
        for i in xrange(-int(window_size/float(2)),int(window_size/float(2))+1):
            if i == 0:
                continue

            actual_aa = self.no + i

            for aa in aa_order:
                if i > 0:
                    name = "evo_%s+%d" %(aa, i)
                else:
                    name = "evo_%s%d" % (aa, i)
                attribute_names.append(name)
                if actual_aa <= 0 or actual_aa > len(self.parent.sequence):
                    attribute_vals.extend([0]*20)
                else:
                    attribute_vals.extend(self.parent.residues[actual_aa].features['PSSM'])

        for idx, val in enumerate(attribute_names):
            self.pro_features[val] = attribute_vals[idx]

        return attribute_names

    def set_termini(self, cutoff=0.1):
        '''
        Sets a binary feature on whether central lysine residue is located in N- or C- terminus
        :param cutoff: Cutoff determining the N- or C- terminus
        :return: Name of the calculated feature
        '''
        if self.no <= math.floor(len(self.parent.sequence)*cutoff) or \
                        self.no >= math.floor(len(self.parent.sequence)-len(self.parent.sequence)*cutoff):
            self.pro_features['termini'] = 1
        else:
            self.pro_features['termini'] = 0
        return ['termini']

    def set_iupred_real(self, window_size):
        '''
        Sets IUPred disorder prediction as the mean disorder predictions of the window
        :param window_size: Size of the window
        :return: Name of the calculated features
        '''
        disorder_real = []
        real_value = 0
        needed_list = []
        window_left = self.no-window_size/2

        #print window_left
        while window_left <= 0:
            needed_list.append(0)
            window_left += 1
            left = 1
        #print f_prob
        left = (window_size/2) - len(needed_list)
        #print left
        for i in range(left):
            #left side
            disorder_real_window = self.parent.residues[self.no-left+i].features['iupred']
            disorder_real.append(disorder_real_window)
            real_value += disorder_real_window
        #residue
        disorder_real_window = self.features['iupred']
        disorder_real.append(disorder_real_window)
        real_value += disorder_real_window

        #right side
        window_right = self.parent.sequence[self.no:self.no+(window_size/2)]
        window_right_len = len(window_right)
        while window_right_len > 0:
            for i in range(window_right_len):
                disorder_real_window = self.parent.residues[self.no+i].features['iupred']
                disorder_real.append(disorder_real_window)
                real_value += disorder_real_window
                window_right_len -= 1
        mean_value = (real_value/len(disorder_real))

        self.pro_features['Window-Disorder-Real'] = mean_value
        return ['Window-Disorder-Real']

    def set_iupred_binary(self, window_size, cutoff=0.5):
        '''
        Sets IUPred disorder prediction as a binary value with a cutoff set to mean disorder prediction of window
        :param window_size: Size of the window
        :param cutoff: Cutoff determining whether window is disordered or not
        :return: Name of the calculated features
        '''
        disorder_real = []
        real_value = 0

        needed_list = []
        window_left = self.no-window_size/2

        while window_left <= 0:
            needed_list.append(0)
            window_left += 1

        #left side
        left = (window_size/2) - len(needed_list)
        for i in range(left):
            disorder_real_window = self.parent.residues[self.no-left+i].features['iupred']
            disorder_real.append(disorder_real_window)
            real_value += disorder_real_window

        #right side
        window_right = self.parent.sequence[self.no:self.no+(window_size/2)]
        window_right_len = len(window_right)
        while window_right_len > 0:
            for i in range(window_right_len):
                disorder_real_window = self.parent.residues[self.no+i].features['iupred']
                disorder_real.append(disorder_real_window)
                real_value += disorder_real_window
                window_right_len -= 1
        mean_value = (real_value/len(disorder_real))
        if mean_value >= cutoff:
            self.pro_features['Window-Disorder-Binary'] = 1
        else:
            self.pro_features['Window-Disorder-Binary'] = 0

        return ['Window-Disorder-Binary']

    def set_hydrophobicity(self, window_size):
        '''
        Sets hydrophobicity features
        :param window_size:
        :return: Names of calculated features
        '''
        attribute_names = []
        attribute_vals = []
        negative_attribute_vals = []
        positive_attribute_vals = []

        total_negative_vals = 0
        total_positive_vals = 0
        total_attribute_vals = 0

        window_start_ind = self.no-int(window_size/2)
        window_start_counter = 1

        window_end_ind = window_start_ind+window_size-1
        window_end_counter = 1
        if window_start_ind<1:
            needed = int(window_size/2)-self.no+1
            window_start_ind += int((math.fabs(window_start_ind)+1))

            for i in xrange(1,int(needed)+1):
                attribute_names.append('w-%s_Hopp-Hydro'%((-int(window_size/float(2))-1)+window_start_counter))
                negative_attribute_vals.append(3.00)
                attribute_vals.append(3.00)
                total_attribute_vals += 3.00
                total_negative_vals += 3.00
                window_start_counter += 1

        if window_end_ind>self.parent.residues[-1].no:
            window_end_ind -= (window_end_ind-self.parent.residues[-1].no)

        for pos in xrange(window_start_ind,window_end_ind+1):

            if pos == self.no:
                continue

            if pos<self.no:
                attribute_names.append('w-%s_Hopp-Hydro'%((-int(window_size/float(2))-1)+window_start_counter))
                neg_attr_val = hopp_hydrophobicity[self.parent.residues[pos].aa]
                negative_attribute_vals.append(neg_attr_val)
                attr_val = hopp_hydrophobicity[self.parent.residues[pos].aa]
                attribute_vals.append(attr_val)
                total_attribute_vals += attr_val
                total_negative_vals += neg_attr_val
                window_start_counter += 1

            elif pos>self.no:
                attribute_names.append('w+%s_Hopp-Hydro'%(window_end_counter))
                pos_attr_val = hopp_hydrophobicity[self.parent.residues[pos].aa]
                positive_attribute_vals.append(pos_attr_val)
                attr_val = hopp_hydrophobicity[self.parent.residues[pos].aa]
                attribute_vals.append(attr_val)
                total_attribute_vals += attr_val
                total_positive_vals += pos_attr_val

                window_end_counter += 1

        if window_end_counter<int(window_size/2)+1:
            needed = int(window_size/2)-window_end_counter+1
            for i in xrange(1,needed+1):
                attribute_names.append('w+%s_Hopp-Hydro' % (window_end_counter))
                positive_attribute_vals.append(3.00)
                attribute_vals.append(3.00)
                total_positive_vals += 3.00
                total_attribute_vals += 3.00

                window_end_counter += 1

        for idx, val in enumerate(attribute_names):
            self.pro_features[val] = attribute_vals[idx]

        avg_attribute_values = total_attribute_vals/len(attribute_vals)

        if len(negative_attribute_vals) > 0:
            avg_negative_attribute_vals = (total_negative_vals/len(negative_attribute_vals))
        else:
            avg_negative_attribute_vals = 0

        if len(positive_attribute_vals) > 0:
            avg_positive_attribute_vals = (total_positive_vals/len(positive_attribute_vals))
        else:
            avg_positive_attribute_vals = 0

        self.pro_features['Negative-Window-Hydro-Avg'] = avg_negative_attribute_vals
        self.pro_features['Positive-Window-Hydro-Avg'] = avg_positive_attribute_vals
        self.pro_features['Window-Hydro-Avg'] = avg_attribute_values
        attribute_names.extend(['Negative-Window-Hydro-Avg', 'Positive-Window-Hydro-Avg', 'Window-Hydro-Avg'])

        return attribute_names

    def set_compositions(self, grouping, window_size, k=1):
        '''
        Sets amino acid composition features in k-mers (e.g. k=1 A, N, k=2 AA, AN, AC,...)
        :param grouping: Amino acid groping (no for no grouping, 'sezerman' for sezerman grouping
        :param window_size: Size of the window
        :param k: Size of k-mers
        :return: Names of the calculated features
        '''
        rawamino_acids = ['A', 'N', 'C', 'Q', 'H', 'L', 'M', 'P', 'T', 'Y', 'R', 'D', 'E', 'G', 'I', 'K', 'F', 'S', 'W', 'V']
        rawamino_acids_sz = ['[IVLM]', '[RKH]', '[DE]', '[QN]', '[ST]', '[A]', '[G]', '[W]', '[C]', '[YF]', '[P]']
        attribute_names = []
        attribute_vals = []
        binary_no = [0]*int(math.pow(20,k))
        binary_sz = [0]*int(math.pow(len(rawamino_acids_sz),k))
        amino_acids = [''.join(i) for i in itertools.product(rawamino_acids, repeat=k)]
        amino_acids_sz = [''.join(i) for i in itertools.product(rawamino_acids_sz, repeat=k)]

        sz_pattern = re.compile(r'\[\w+?\]')

        window_start_ind = self.no-int(window_size/2)
        window_start_counter = 1

        window_end_ind = window_start_ind+window_size-1
        window_end_counter = 1
        if window_start_ind<1:
            needed = int(window_size/2)-self.no+1
            window_start_ind += int((math.fabs(window_start_ind)+1))

            for i in xrange(1,int(needed)+1):
                if grouping=="no":
                    for aa in amino_acids:
                        attribute_names.append('w-%s_%s'%(aa, (-int(window_size/float(2))-1)+window_start_counter))
                    attribute_vals.extend(binary_no)
                elif grouping=="sezerman":
                    for aa in amino_acids_sz:
                        attribute_names.append('w-%s_%s-SZ'%(aa, (-int(window_size/float(2))-1)+window_start_counter))
                    attribute_vals.extend(binary_sz)
                window_start_counter += 1

        if window_end_ind>self.parent.residues[-1].no:
            window_end_ind -= (window_end_ind-self.parent.residues[-1].no)

        for pos in xrange(window_start_ind,window_end_ind+1):
            if pos == self.no:
                continue
            if pos<(self.no-(k-1)):
                if grouping=="no":
                    temp = list(binary_no)
                    tmpstr = ""
                    for i in xrange(0,k):
                        tmpstr += self.parent.residues[pos+i].aa
                    indexofstr = amino_acids.index(tmpstr)
                    temp[indexofstr] = 1

                    for aa in amino_acids:
                        attribute_names.append('w-%s_%s'%(aa, (-int(window_size/float(2))-1)+window_start_counter))
                    attribute_vals.extend(temp)

                elif grouping == "sezerman":
                    temp = list(binary_sz)
                    tmpstr = ""
                    for i in xrange(0,k):
                        tmpstr += "["+sz_grouping[self.parent.residues[pos+i].aa]+"]"
                    indexofstr = amino_acids_sz.index(tmpstr)
                    temp[indexofstr] = 1

                    for aa in amino_acids_sz:
                        attribute_names.append('w-%s_%s-SZ'%(aa, (-int(window_size/float(2))-1)+window_start_counter))
                    attribute_vals.extend(temp)

                window_start_counter += 1
            elif pos > self.no and pos <= (self.no+int(window_size/2)-(k-1)):
                if grouping=="no":
                    temp = list(binary_no)
                    tmpstr = ""
                    for i in xrange(0,k):
                        if (pos+i) < len(self.parent.residues):
                            tmpstr += self.parent.residues[pos+i].aa
                        else:
                            break
                    if len(tmpstr) < k:
                        break
                    indexofstr = amino_acids.index(tmpstr)
                    temp[indexofstr] = 1

                    for aa in amino_acids:
                        attribute_names.append('w+%s_%s'%(aa, window_end_counter))
                    attribute_vals.extend(temp)

                elif grouping=="sezerman":
                    temp = list(binary_sz)
                    tmpstr = ""
                    for i in xrange(0,k):
                        if (pos+i) < len(self.parent.residues):
                            tmpstr += "["+sz_grouping[self.parent.residues[pos+i].aa]+"]"
                        else:
                            break
                    items = sz_pattern.findall(tmpstr)
                    if len(items) < k:
                        break
                    indexofstr = amino_acids_sz.index(tmpstr)
                    temp[indexofstr] = 1

                    for aa in amino_acids_sz:
                        attribute_names.append('w+%s_%s-SZ'%(aa, window_end_counter))
                    attribute_vals.extend(temp)

                window_end_counter += 1

        if (window_end_counter-1)<(int(window_size/2)-(k-1)):
            needed = int(window_size/2)-(k-1)-(window_end_counter-1)
            for i in xrange(1,needed+1):
                if grouping=="no":
                    for aa in amino_acids:
                        attribute_names.append('w+%s_%s'%(aa, window_end_counter))
                    attribute_vals.extend(binary_no)
                elif grouping=="sezerman":
                    for aa in amino_acids_sz:
                        attribute_names.append('w+%s_%s-SZ'%(aa, window_end_counter))
                    attribute_vals.extend(binary_sz)
                window_end_counter += 1

        for idx, val in enumerate(attribute_names):
            self.pro_features[val] = attribute_vals[idx]
        return attribute_names

    def get_window(self, window_size):
        '''
        Returns window sequence
        :param window_size: Size of the sequence window
        :return: Tuple containing left and right side sequences of the central lysine residue
        '''
        resloc = self.no-1
        left = ""
        if (resloc-((window_size-1)/2) < 0):
            missing = int(math.fabs(resloc-((window_size-1)/2)))
            left += '-'*missing
            left += self.parent.sequence[0:resloc]
        else:
            left = self.parent.sequence[(resloc-((window_size-1)/2)):resloc]

        right = ""
        if (resloc+1+((window_size-1)/2) > len(self.parent.sequence)):
            missing = (resloc+((window_size-1)/2)+1) - len(self.parent.sequence)
            right = self.parent.sequence[resloc+1:]
            right += '-'*missing
        else:
            right = self.parent.sequence[resloc+1:(resloc+1+((window_size-1)/2))]

        return (left, right)

    def get_instance(self, feats, scaler):
        '''
        Wrapper for encoding the window for SVM prediction
        :param feats: features to be used in the encoding
        :param scaler: classifier scaler object
        :return: Bunch object
        '''
        x = []
        y = 0

        for attr in feats:
            if not attr in self.pro_features.keys() or self.pro_features[attr] == '?':
                x.append(numpy.nan)
            else:
                x.append(self.pro_features[attr])

        if self.features["class"]:
            y = [1]
        target_names = ['0', '1']

        return Bunch(data=scaler.transform(x), feature_names=feats, target=y, target_names=target_names)