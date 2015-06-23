from __future__ import absolute_import
import numpy
from sklearn.metrics import *
from sklearn.preprocessing import Imputer
import re
import uuid
import sys
import copy
import math
import logging
import itertools
import exceptions as exceptions

class Bunch(dict):
    """Container object for datasets: dictionary-like object that
       exposes its keys as attributes."""

    def __init__(self, **kwargs):
        dict.__init__(self, kwargs)
        self.__dict__ = self

psipred_dict = { 'H' : 0, 'E' : 1, 'C': 2}

hopp_hydrophobicity = { 'A' : -0.50, 'C' : -1.00, 'D' : 3.00, 'E' : 3.00, 'F' : -2.50, 'G' : 0.00, 'H' : -0.50,
                        'I' : -1.80, 'K' : 3.00, 'L' : -1.80, 'M' : -1.30, 'N' : 0.20, 'P' : 0.00, 'Q' : 0.20,
                        'R' : 3.00,  'S' : 0.30, 'T' : -0.40, 'V' : -1.50, 'W' : -3.40, 'Y' : -2.30 }

aa_1to3 = { "A" : "ALA", "R" : "ARG", "N" : "ASN", "D" : "ASP", "C" : "CYS", "E" : "GLU", "Q" : "GLN", "G" : "GLY",
            "H" : "HIS", "I" : "ILE", "L" : "LEU", "K" : "LYS", "M" : "MET", "F" : "PHE", "P" : "PRO", "S" : "SER",
            "T" : "THR", "W" : "TRP", "Y" : "TYR", "V" : "VAL" }

aa_3to1 = { "ALA" : "A", "ARG" : "R", "ASN" : "N", "ASP" : "D", "CYS" : "C", "GLU" : "E", "GLN" : "Q", "GLY" : "G",
            "HIS" : "H", "ILE" : "I", "LEU" : "L", "LYS" : "K", "MET" : "M", "PHE" : "F", "PRO" : "P", "SER" : "S",
            "THR" : "T", "TRP" : "W", "TYR" : "Y", "VAL" : "V" }

# Normal Amino Acid Alphabet
aa_dict = { "ALA" : 0, "ARG" : 1, "ASN" : 2, "ASP" : 3, "CYS" : 4, "GLN" : 5, "GLU" : 6, "GLY" : 7, "HIS" : 8,
            "ILE" : 9, "LEU" : 10, "LYS" : 11, "MET" : 12, "PHE" : 13, "PRO" : 14, "SER" : 15, "THR" : 16, "TRP" : 17,
            "TYR" : 18, "VAL" : 19}

# Sezerman Grouping
sz_aa_dict = { "ALA" : 5, "ARG" : 1, "ASN" : 3, "ASP" : 2, "CYS" : 8, "GLN" : 3, "GLU" : 2, "GLY" : 6, "HIS" : 1,
               "ILE" : 0, "LEU" : 0, "LYS" : 1, "MET" : 0, "PHE" : 9, "PRO" : 10, "SER" : 4, "THR" : 4, "TRP" : 7,
               "TYR" : 9, "VAL" : 0}

sz_grouping = { "A" : "A", "R" : "RKH", "N" : "QN", "D" : "DE", "C" : "C", "E" : "DE", "Q" : "QN", "G" : "G",
            "H" : "RKH", "I" : "IVLM", "L" : "IVLM", "K" : "RKH", "M" : "IVLM", "F" : "YF", "P" : "P", "S" : "ST",
            "T" : "ST", "W" : "W", "Y" : "YF", "V" : "IVLM"}

volumes = { 'G' : 37.5, 'A' : 54.7, 'V' : 85.05, 'L' : 101.92, 'I' : 99.86, 'P' : 74.8, 'F' : 116.06, 'Y' : 118,
            'W' : 138.16, 'H' : 93.06, 'M' : 99.62, 'C' : 67.71, 'S' : 54.95, 'T' : 71.12, 'N' : 71.56, 'Q' : 88.16,
            'D' : 68.56, 'E' : 84, 'K' : 105.7, 'R' : 121.2 }
class Residue:
    # Residue wrapper for protein sequences or structures.
    parent = None   # Parent sequence
    no = 0   # Amino acid no in sequence
    aa = ''   # Amino acid
    features = {}   # Unprocessed feature dictionary
    pro_features = None   # Processed feature pointer

    def __init__(self, no=0, aa='', parent=None, raw_features=None, **kwargs):
        self.parent = parent
        self.no = no
        self.aa = aa
        if raw_features:
            self.features = raw_features
        else:
            self.features = {}
        self.pro_features = {}


    def set_all_seq(self, window_size=17):
        flex_attr = self.set_flexibility(window_size=window_size)
        comp1_attr = self.set_compositions(grouping="no", k=1, window_size=window_size)
        comp1_attr_sz = self.set_compositions(grouping="sezerman", k=1, window_size=window_size)
        iupr_attr = self.set_iupred_real(window_size=window_size)
        iupb_attr = self.set_iupred_binary(window_size=window_size)
        psip_attr = self.set_secondarystructure(sstype='psipred')
        psipss_attr = self.set_ssends(source='psipred')
        sqend_attr = self.set_seqends()
        termini_attr = self.set_termini()
        winfo_attr = self.set_window_info(grouping="hopp", window_size=window_size)
        aacounts_attr = self.set_aacounts(grouping="no", window_size=window_size)
        aacounts_sz_attr = self.set_aacounts(grouping="sezerman", window_size=window_size)
        aafreqs_attr = self.set_aafreqratios(grouping="no", window_size=window_size)
        aafreqs_sz_attr = self.set_aafreqratios(grouping="sezerman", window_size=window_size)
        #special_attr = self.set_specials(window_size=window_size)
        pssm_attr = self.set_pssm(window_size=window_size)

        return {'FLEX':flex_attr, 'COMP1': comp1_attr,
                 'COMP1_SZ': comp1_attr_sz, 'PSSM': pssm_attr,
                 'AACNT': aacounts_attr, 'AACNTSZ': aacounts_sz_attr,
                 'AAFRQ': aafreqs_attr, 'AAFRQSZ': aafreqs_sz_attr,
                 'IUPR': iupr_attr, 'IUPB': iupb_attr, 'PSIP': psip_attr,
                 'PSIPSS': psipss_attr, 'SQEND': sqend_attr, 'TERMINI': termini_attr,
                 'WINFO': winfo_attr #'SPECIAL': special_attr
        }

    def set_selected(self,window_size=21):
        comp1_attr = self.set_compositions(grouping="no", k=1, window_size=window_size)
        comp1_attr_sz = self.set_compositions(grouping="sezerman", k=1, window_size=window_size)
        iupr_attr = self.set_iupred_real(window_size=window_size)
        iupb_attr = self.set_iupred_binary(window_size=window_size)
        sqend_attr = self.set_seqends()
        termini_attr = self.set_termini()
        winfo_attr = self.set_window_info(grouping="hopp", window_size=window_size)
        aacounts_attr = self.set_aacounts(grouping="no", window_size=window_size)
        aacounts_sz_attr = self.set_aacounts(grouping="sezerman", window_size=window_size)
        aafreqs_attr = self.set_aafreqratios(grouping="no", window_size=window_size)
        aafreqs_sz_attr = self.set_aafreqratios(grouping="sezerman", window_size=window_size)
        #special_attr = self.set_specials(window_size=window_size)
        pssm_attr = self.set_pssm(window_size=window_size)

        return { 'COMP1': comp1_attr,
                 'COMP1_SZ': comp1_attr_sz, 'PSSM': pssm_attr,
                 'AACNT': aacounts_attr, 'AACNTSZ': aacounts_sz_attr,
                 'AAFRQ': aafreqs_attr, 'AAFRQSZ': aafreqs_sz_attr,
                 'IUPR': iupr_attr, 'IUPB': iupb_attr,
                 'SQEND': sqend_attr, 'TERMINI': termini_attr,
                 'WINFO': winfo_attr #'SPECIAL': special_attr
        }

    def set_aacounts(self, grouping="no", window_size=7):
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
        else:
            raise exceptions.InvalidGrouping("Invalid Amino Acid Grouping: %s" % grouping)

        for idx, val in enumerate(attribute_names):
            self.pro_features[val] = attribute_vals[idx]

        return attribute_names

    def set_aafreqratios(self, grouping="no", window_size=7):
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
        else:
            raise exceptions.InvalidGrouping("Invalid Amino Acid Grouping: %s" % grouping)

        for idx, val in enumerate(attribute_names):
            self.pro_features[val] = attribute_vals[idx]

        return attribute_names

    def set_pssm(self, window_size=7):
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
        if self.no <= math.floor(len(self.parent.sequence)*0.1) or self.no >= math.floor(len(self.parent.sequence)-len(self.parent.sequence)*0.1):
            self.pro_features['termini'] = 1
        else:
            self.pro_features['termini'] = 0
        return ['termini']

    def set_seqends(self, cutoff=0.1):
        if self.no <= math.floor(len(self.parent.sequence)*0.1) or self.no >= math.floor(len(self.parent.sequence)-len(self.parent.sequence)*0.1):
            self.pro_features['seqlen'] = (-1)*len(self.parent.sequence)
        else:
            self.pro_features['seqlen'] = len(self.parent.sequence)
        return ['seqlen']


    def set_flexibility(self, window_size):
        f_prob = []
        r_prob = []
        needed_list = []
        total_p = 0
        total_r = 0
        window_left = self.no-window_size/2

        while window_left <= 0:
            needed = 0 - window_left + 1

            needed_list.append(0)
            window_left += 1
            left = 1

        left = (window_size/2) - len(needed_list)


        for i in range(left):

            if self.parent.features['flexpred'][self.no-left+i] == "F":
                flexible = self.parent.residues[self.no-left+i].features['flexpredconf']
                f_prob.append(flexible)
                total_p += float(flexible)
            elif self.parent.features['flexpred'][self.no-left+i] == "R":
                rigidity = self.parent.residues[self.no-left+i].features['flexpredconf']
                r_prob.append(rigidity)
                total_r += float(rigidity)

        if self.features['flexpred'] == 'F':
            flexible = self.features['flexpredconf']
            f_prob.append(flexible)
            total_p += float(flexible)
        else:
            rigidity = self.features['flexpredconf']
            r_prob.append(rigidity)
            total_r += float(rigidity)

        #right side
        window_right = self.parent.sequence[self.no:self.no+(window_size/2)]
        window_right_len = len(window_right)
        while window_right_len > 0:
            for i in range(window_right_len):
                if self.parent.features['flexpred'][self.no+i] == "F":
                    flexible = self.parent.residues[self.no+i].features['flexpredconf']
                    f_prob.append(flexible)
                    total_p += float(flexible)
                elif self.parent.features['flexpred'][self.no+i] == "R":
                    rigidity = self.parent.residues[self.no+i].features['flexpredconf']
                    r_prob.append(rigidity)
                    total_r += float(rigidity)
                window_right_len -= 1

        ratio = 0
        if total_r == 0:
            ratio = 1
        else:
            if total_p / total_r > 1:
                ratio = 1
            else:
                ratio = 0
        self.pro_features['Window-Flexibility-ConfRat'] = ratio
        if len(f_prob) >= (window_size/2):
            self.pro_features['Window-Flexibility'] = 1
        else:
            self.pro_features['Window-Flexibility'] = 0
        return ['Window-Flexibility-ConfRat', 'Window-Flexibility']

    # def set_specials(self, window_size):
    #     window = ''
    #     if (self.no-(window_size/2)+1) < 0:
    #         missing = int(math.fabs(self.no-(window_size/2)+1))
    #         window += 'X'*missing
    #         window += self.parent.sequence[0:(self.no-1)]
    #
    #     else:
    #         window += self.parent.sequence[(self.no-(window_size/2)+1):(self.no-1)]
    #
    #     if (self.no+((window_size-1)/2)) > len(self.parent):
    #         needed = int(math.fabs(self.no+((window_size-1)/2)-len(self.parent)))
    #         window += self.parent.sequence[(self.no):]
    #         window += 'X'*needed
    #     else:
    #         window += self.parent.sequence[(self.no):self.no+((window_size-1)/2)]
    #
    #     if 'D' in window or 'E' in window:
    #         self.pro_features['wDE'] = 1
    #     else:
    #         self.pro_features['wDE'] = 0
    #
    #     if 'K' in window:
    #         self.pro_features['wK'] = 1
    #     else:
    #         self.pro_features['wK'] = 0
    #         return ['wDE', 'wK']

    def set_secondarystructure(self,sstype):
        attribute_names = []
        attribute_vals = []


        attribute_names = ['pHelix', 'pBeta', 'pCoil']
        attribute_vals = [0]*3
        attribute_vals[psipred_dict[self.features['psipred']]] = 1

        for idx, val in enumerate(attribute_names):
            self.pro_features[val] = attribute_vals[idx]
        return ['pHelix', 'pBeta', 'pCoil']

    def set_ssends(self, source="dssp", cutoff=0.15):
        attribute_names = ['AtSSEnd_'+source]
        attribute_vals = [0]
        helices = self.parent.features["helices"+source]
        strands = self.parent.features["strands"+source]

        def is_in_ends(ends, position, cutoff = 0.15):
            in_end = False
            for sstr in ends:
                length = sstr[1]-sstr[0]+1
                cut_line = round(length*cutoff)
                if ((position<(sstr[0]+cut_line) or position>(sstr[1]-cut_line)) and (position<=sstr[1] and position>=sstr[0])):
                    in_end = True
            if (in_end):
                return 1
            else:
                return 0

        ss_val = psipred_dict[self.features["psipred"]]

        if (ss_val==1):
            attribute_vals[0] = is_in_ends(strands, self.no, cutoff)
        elif (ss_val==0):
            attribute_vals[0] = is_in_ends(helices, self.no, cutoff)

        for idx, val in enumerate(attribute_names):
            self.pro_features[val] = attribute_vals[idx]
        return ['AtSSEnd_psipred']

    def set_wesa(self):
        if 'wesa' in self.features.keys():
            self.pro_features['SolAcc-WESA'] = self.features["wesa"]
        else:
            raise exceptions.NoFeature('WESA features are not available for this residue.')

        return ['SolAcc-WESA']

    def set_iupred_real(self, window_size):
        disorder_real = []
        real_value = 0
        needed_list = []
        window_left = self.no-window_size/2
        #print window_left
        while window_left <= 0:
            needed = 0 - window_left + 1
            #print "Needed: %s" % needed
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

    def set_iupred_binary(self, window_size):

        disorder_real = []
        real_value = 0
        #positive residues
        #if 'iupred' in self.features.keys():
        #    self.pro_features['DisorderReal'] = self.features["iupred"]
        #else:
        #    raise NoFeature('IUPRED features are not available for this residue.')
        needed_list = []
        window_left = self.no-window_size/2
        #print window_left
        while window_left <= 0:
            needed = 0 - window_left + 1
            #print "Needed: %s" % needed
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
        if mean_value >= 0.5:
            self.pro_features['Window-Disorder-Binary'] = 1
        else:
            self.pro_features['Window-Disorder-Binary'] = 0

        #if 'iupred' in self.features.keys():
        #    if self.features["iupred"] > 0.5:
        #        self.pro_features['DisorderBinary'] = 1
        #    else:
        #        self.pro_features['DisorderBinary'] = 0
        #else:
        #    raise NoFeature('IUPRED features are not available for this residue.')

        return ['Window-Disorder-Binary']


    def set_window_info(self, grouping, window_size):

            attribute_names = []
            negative_attribute_names = []
            positive_attribute_names = []
            attribute_vals = []
            negative_attribute_vals = []
            positive_attribute_vals = []

            total_negative_vals = 0
            total_positive_vals = 0
            total_attribute_vals = 0

            before_vol = 0.
            after_vol = 0.

            window_start_ind = self.no-int(window_size/2)
            window_start_counter = 1

            window_end_ind = window_start_ind+window_size-1
            window_end_counter = 1
            if window_start_ind<1:
                needed = int(window_size/2)-self.no+1
                window_start_ind += int((math.fabs(window_start_ind)+1))

                for i in xrange(1,int(needed)+1):
                    if grouping=="hopp":
                        #negative_attribute_names.append('w-%s_Hopp-Hydro'%(window_start_counter))
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
                    #attribute_names.extend(['w-%s_AA'%(window_start_counter), 'w-%s_Hydro'%(window_start_counter)])
                    #attribute_vals.extend([seq[pos].aa_name, hopp_hydrophobicity[seq[pos].aa_name]])
                    if grouping=="hopp":
                        #negative_attribute_names.append('w-%s_Hopp-Hydro'%(window_start_counter))
                        attribute_names.append('w-%s_Hopp-Hydro'%((-int(window_size/float(2))-1)+window_start_counter))
                        neg_attr_val = hopp_hydrophobicity[self.parent.residues[pos].aa]
                        negative_attribute_vals.append(neg_attr_val)
                        attr_val = hopp_hydrophobicity[self.parent.residues[pos].aa]
                        attribute_vals.append(attr_val)
                        total_attribute_vals += attr_val
                        total_negative_vals += neg_attr_val

                    before_vol += volumes[self.parent.residues[pos].aa]
                    window_start_counter += 1

                elif pos>self.no:
                    if grouping=="hopp":
                        #positive_attribute_names.append('w+%s_Hopp-Hydro'%(window_end_counter))
                        attribute_names.append('w+%s_Hopp-Hydro'%(window_end_counter))
                        pos_attr_val = hopp_hydrophobicity[self.parent.residues[pos].aa]
                        positive_attribute_vals.append(pos_attr_val)
                        attr_val = hopp_hydrophobicity[self.parent.residues[pos].aa]
                        attribute_vals.append(attr_val)
                        total_attribute_vals += attr_val
                        total_positive_vals += pos_attr_val

                    after_vol += volumes[self.parent.residues[pos].aa]
                    window_end_counter += 1

            if window_end_counter<int(window_size/2)+1:
                needed = int(window_size/2)-window_end_counter+1
                for i in xrange(1,needed+1):
                    if grouping=="hopp":
                        #positive_attribute_names.append('w+%s_Hopp-Hydro' % (window_end_counter))
                        attribute_names.append('w+%s_Hopp-Hydro' % (window_end_counter))
                        positive_attribute_vals.append(3.00)
                        attribute_vals.append(3.00)
                        total_positive_vals += 3.00
                        total_attribute_vals += 3.00

                    window_end_counter += 1

            attribute_names.extend(['BeforeVol', 'AfterVol', 'Difference'])
            attribute_vals.extend([before_vol, after_vol, (after_vol-before_vol)])
            #if len(attribute_vals) != 9 or not ('w+3_Hopp-Hydro' in attribute_names):
            #    print attribute_names
            for idx, val in enumerate(attribute_names):
                self.pro_features[val] = attribute_vals[idx]
            for idx, val in enumerate(negative_attribute_names):
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

                elif grouping=="sezerman":
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
            elif pos>self.no and pos<=(self.no+int(window_size/2)-(k-1)):
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

    # def get_instance(self, feats, imputer, scaler):
    #     x = []
    #     y = 0
    #
    #     for attr in feats:
    #         if not attr in self.pro_features.keys() or self.pro_features[attr] == '?':
    #             x.append(numpy.nan)
    #         else:
    #             x.append(self.pro_features[attr])
    #
    #     if self.features["class"]:
    #         y = [1]
    #     target_names = ['0', '1']
    #
    #     return Bunch(data=scaler.transform(imputer.transform(x)), feature_names=feats, target=y, target_names=target_names)
    # def get_instance(self, encoding, tags, form='scikit'):
    #     if form == 'svm':
    #         x = {}
    #         y = []
    #         attr_list = encoding.get_attrs(tags, form='svm') # w/o class!
    #         for i, attr in enumerate(attr_list):
    #             x[i] = self.pro_features[attr]
    #         if self.features["class"]:
    #             return (1, x)
    #         else:
    #             return (0, x)
    #     elif form == 'csv':
    #         x = ""
    #         y = []
    #         attr_list = encoding.get_attrs(tags, form='csv') # w/o class!
    #         for i, attr in enumerate(attr_list):
    #             x += str(self.pro_features[attr]) + ","
    #         if self.features["class"]:
    #             return (1, x[:-1])
    #         else:
    #             return (0, x[:-1])
    #     elif form == 'scikit':
    #         x = []
    #         y = 0
    #         for tag in tags:
    #             for attr in encoding.sets[tag]:
    #                 if not attr in self.pro_features.keys() or self.pro_features[attr] == '?':
    #                     x.append(numpy.nan)
    #                 else:
    #                     x.append(self.pro_features[attr])
    #         if self.features["class"]:
    #             y = [1]
    #         (desc, feature_names) = encoding.get_features(tags=tags, form='scikit')
    #         target_names = ['0', '1']
    #
    #         return Bunch(data=[x], feature_names=feature_names, DESC=desc, target=y, target_names=target_names)