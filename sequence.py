from __future__ import print_function
from __future__ import absolute_import

import re
import os
import sys
import tempfile
import logging
import logging.config
import random
import numpy
import residue as residue
import exceptions as exceptions


#PSIPRED Prediction alphabet
psipred_dict = { 'H' : 0, 'E' : 1, 'C': 2}

# Normal Amino Acid Alphabet
aa_dict = {"ALA": 0, "ARG": 1, "ASN": 2, "ASP": 3, "CYS": 4, "GLN": 5, "GLU": 6, "GLY": 7, "HIS": 8,
           "ILE": 9, "LEU": 10, "LYS": 11, "MET": 12, "PHE": 13, "PRO": 14, "SER": 15, "THR": 16, "TRP": 17,
           "TYR": 18, "VAL": 19}

class Seq(object):
    '''
        Sequence and structure wrapper.
    '''
    identifier = ''
    sequence = ''
    residues = []
    features = {}
    organism = None

    def __init__(self, identifier, sequence, features=None, verbosity=0, **kwargs):
        self.identifier = identifier
        self.sequence = sequence
        self.residues = []
        self.residues.append(None)
        self.verbosity = 0
        if 'organism' in kwargs:
            self.organism = kwargs['organism']
        for i in xrange(0, len(self.sequence)):
            tmp = residue.Residue(no=(i + 1), aa=self.sequence[i], parent=self)
            tmp.features["class"] = False
            tmp.features["hbonds"] = []
            self.residues.append(tmp)
        if features:
            self.features = features
        else:
            self.features = {}



    def add_pssm(self, filename=''):
            self.parse_from_pssm(filename)

    def add_flexpred(self, filename):
        self.parse_from_flexpred(filename)

    def add_psipred(self, filename):
        self.parse_from_psipred(filename)

    def add_iupred(self, filename=''):
            self.parse_from_iupred(filename)

    def add_wesa(self, filename):
        self.parse_from_wesa(filename)

    def add_sites(self, sites = None):
        self.parse_sites_from_list(sites)

    def parse_from_iupred(self, filename):
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
                    raise exceptions.InvalidFile('Mismatch in IUPRED File with Sequence!')

    def parse_sites_from_list(self, sites):
        if not "sites" in self.features.keys():
            self.features["sites"] = []

        for site in sites:
            self.features["sites"].append(site)
            self.residues[site].features["class"] = True

    def parse_from_pssm(self, ifile):
        ifh = open(ifile, "r")
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
                        raise exceptions.InvalidFile("Invalid PSSM File! Amino acid mismatch at line %d for aminoacid %d in file %s.\nHave "
                                      "%s, got %s." % (idx, rc, ifile, self.residues[rc].aa, aa))

                    elems = line_match.group(3).split()[0:20]
                    self.residues[rc].features['PSSM'] = [ int(x) for x in elems ]
                else:
                    continue

        self.features['PSSM_AA_Order'] = aa_order

    def parse_from_psipred(self, ofile):
        ifh = open(ofile, "r")
        lines = ifh.readlines()
        ifh.close()

        confre = re.compile(r'Conf:\s*(.*)')
        predre = re.compile(r'Pred:\s*(.*)')

        confseq = ""
        predseq = ""

        for line in lines:
            cm = confre.search(line)
            pm = predre.search(line)
            if (cm):
                confseq += cm.group(1)
            elif (pm):
                predseq += pm.group(1)

        for pos in xrange(1, len(predseq) + 1):
            self.residues[pos].features["psipred"] = predseq[pos - 1]
            self.residues[pos].features["psipredconf"] = int(confseq[pos - 1])

        self.features["psipred"] = predseq
        self.features["psipredconf"] = confseq
        self.set_helix_ends()
        self.set_strand_ends()

    def parse_from_flexpred(self, ofile):
        ifh = open(ofile, "r")
        lines = ifh.readlines()
        ifh.close()
        resre = re.compile(r'(\d+)\s*(\w+)\s*(\w+)\s*(\d+\.\d+)')
        flexseq = ""
        for idx, line in enumerate(lines):
            lm = resre.search(line)
            if (lm):
                rc = int(lm.group(1))
                aa = lm.group(2)
                fl = lm.group(3)
                fl_conf = lm.group(4)
                if self.residues[rc].aa == aa:
                    flexseq += fl
                    self.residues[rc].features["flexpred"] = fl
                    self.residues[rc].features["flexpredconf"] = fl_conf
                else:
                    raise exceptions.InvalidFile("Invalid FlexPred File! Amino acid mismatch at line %d for aminoacid %d in file %s.\nHave "
                                      "%s, got %s." % (idx, rc, ofile, self.residues[rc].aa, aa))

        self.features["flexpred"] = flexseq

    def parse_from_wesa(self, ofile):
        ifh = open(ofile, "r")
        lines = ifh.readlines()
        ifh.close()

        for idx, line in enumerate(lines):
            if idx >= 16:
                if line.startswith('*'):
                    break
                elems = line.split()
                if self.residues[int(elems[0])].aa == elems[1]:
                    self.residues[int(elems[0])].features["wesa"] = int(elems[7][0:1])
                else:
                    raise exceptions.InvalidFile("Invalid WESA File. Amino acid mismatch at line %d." % (idx))

    def set_aa_freq(self):
        rawamino_acids = ['A', 'N', 'C', 'Q', 'H', 'L', 'M', 'P', 'T', 'Y', 'R', 'D', 'E', 'G', 'I', 'K', 'F', 'S', 'W', 'V']
        rawamino_acids_sz = ['[IVLM]', '[RKH]', '[DE]', '[QN]', '[ST]', '[A]', '[G]', '[W]', '[C]', '[YF]', '[P]']

        attribute_names = ['freqSeq'+i for i in rawamino_acids]
        attribute_vals = [(self.sequence.count(val))/float(len(self.sequence)) for val in rawamino_acids]

        attribute_names.extend(['freqSeq'+i for i in rawamino_acids_sz])
        attribute_vals.extend([len(re.findall(val, self.sequence))/float(len(self.sequence)) for val in rawamino_acids_sz])

        for idx, val in enumerate(attribute_names):
            self.features[val] = attribute_vals[idx]

    def set_helix_ends(self):
        ss_sequence = []
        if self.features["psipred"]:
            for i in xrange(0,len(self.features["psipred"])):
                ss_sequence.append(psipred_dict[self.features["psipred"][i]])
        else:
            raise exceptions.NoFeature('PSIPRED Secondary Structure is not added to sequence!')

        helices = []
        helix_counter = 0
        helix_found = False
        for i in xrange(0, len(ss_sequence)):
            if (ss_sequence[i] == 0 and not helix_found):
                helix_found = True
                helices.append([i+1])
            elif (helix_found and ss_sequence[i] != 0):
                helices[helix_counter].append(i + 1)
                helix_counter += 1
                helix_found = False
        self.features["helicespsipred"] = helices

    def set_strand_ends(self):
        ss_sequence = []

        if self.features["psipred"]:
            for i in xrange(0,len(self.features["psipred"])):
                ss_sequence.append(psipred_dict[self.features["psipred"][i]])
        else:
            raise exceptions.NoFeature('PSIPRED Secondary Structure is not added to sequence!')

        strands = []
        strand_counter = 0
        strand_found = False
        for i in xrange(0, len(ss_sequence)):
            if (ss_sequence[i] == 1 and not strand_found):
                strand_found = True
                strands.append([i+1])
            elif (strand_found and ss_sequence[i] != 1):
                strands[strand_counter].append(i + 1)
                strand_counter += 1
                strand_found = False
        self.features["strandspsipred"] = strands

    def get_fasta(self):
        '''
            Returns FASTA (Pearson) formatted string of desired sequence type.
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

class Encoding:
    '''
        Processed feature wrapper for residues.
    '''
    counts = [0,0]
    classes = ([],[])
    sets = dict()

    def __init__(self, residues = None, sets = None):
        self.counts = [0, 0]
        self.classes = ([],[])
        if residues:
            for residue in residues:
                if hasattr(residue, 'features'):
                    if residue.features["class"]==True:
                        self.classes[0].append(residue)
                        self.counts[0] += 1
                    else:
                        self.classes[1].append(residue)
                        self.counts[1] += 1
                else:
                    if residue.pro_features["class"]==True:
                        self.classes[0].append(residue)
                        self.counts[0] += 1
                    else:
                        self.classes[1].append(residue)
                        self.counts[1] += 1
        else:
            self.residues = []
        if sets:
            self.sets = sets
        else:
            self.sets = {}

    def add_residue(self, res):
        self.residues.append(res)
        if res.features["class"]==True:
            self.classes[0].append(res)
            self.counts[0] += 1
        else:
            self.classes[1].append(res)
            self.counts[1] += 1

    def extend_residues(self, residues):
        self.residues.extend(residues)

        for residue in residues:
            if residue.features["class"]==True:
                self.classes[0].append(residue)
                self.counts[0] += 1
            else:
                self.classes[1].append(residue)
                self.counts[1] += 1

    def sample_dataset(self, ratio):
        # Negative to positive ratio
        neg_amount = ratio*self.counts[0]
        random.shuffle(self.classes[1])
        if neg_amount >= self.counts[1]:
            sampled = self.classes[1]
        else:
            sampled = self.classes[1][1:neg_amount+1]
        sets = (self.classes[0], sampled)

        return sets

    def get_features(self, tags, form='arff', relation=''):
        header = None
        if relation == '':
                relation = 'Modification_'+''.join(tags)[1:20]
        if form == 'arff':
            header = ''
            header = "@relation "+relation+"\n"
            for tag in tags:
                for attr in self.sets[tag]:
                    real = 0
                    binary = 0
                    attrtype = ""
                    for i in xrange(0,5):
                        ind = random.randint(0, len(self.classes[1])-1)
                        if not attr in self.classes[1][ind].pro_features.keys():
                            i -= 1
                            continue
                        if isinstance(self.classes[1][ind].pro_features[attr], int):
                            binary += 1
                        elif isinstance(self.classes[1][ind].pro_features[attr], float):
                            real += 1
                    if real>binary or attr == 'seqlen':
                        attrtype = 'numeric'
                    elif binary>real and real==0:
                        attrtype = '{0,1}'
                    if attr.startswith('evo'):
                        attrtype = 'numeric'

                    header += "@attribute "+attr+" "+attrtype+"\n"
            header += "@attribute status {0, 1}\n"
            header += "@data\n"
        elif form == 'scikit':
            feature_names = []
            for tag in tags:
                for attr in self.sets[tag]:
                    feature_names.append(attr)
            header = (relation, feature_names)
        return header

    def get_attrs(self, tags, form='arff'):
        if form == 'svm' or form == 'csv':
            attrs = []
            for tag in tags:
                for attr in self.sets[tag]:
                    attrs.append(attr)
            return attrs

    def get_encoding(self, ratio=None, tags=None, randomize=True, form='scikit', imputer=None):
        encoding = None
        instances = []
        if ratio == None:
            sample = self.classes
        else:
            sample = self.sample_dataset(ratio)

        if tags == None:
            tags = self.classes[1][1].features.keys()

        instances.extend(sample[0])
        instances.extend(sample[1])
        if randomize:
            random.shuffle(instances)

        if form=='arff':
            encoding = ""
            for inst in instances:
                instencoded = ""
                for tag in tags:
                    for attr in self.sets[tag]:
                        if not attr in inst.pro_features.keys():
                            instencoded += '?, '
                        else:
                            instencoded += str(inst.pro_features[attr])
                            instencoded += ", "
                if inst.features["class"]:
                    instencoded += "1\n"
                else:
                    instencoded += "0\n"
                encoding += instencoded
        elif form == 'scikit':
            encoded_instances = []
            targets = []
            for inst in instances:
                instance = []
                for tag in tags:
                    for attr in self.sets[tag]:
                        if (not attr in inst.pro_features.keys()) or (inst.pro_features[attr] == '?'):
                            print("We have NaN over hereee!")
                            instance.append(numpy.nan)
                        else:
                            instance.append(inst.pro_features[attr])
                    if len(instance) < len(self.sets[tag]):
                        raise BaseException("We have an issue here! Instance length is not equal to tag length.")
                encoded_instances.append(instance)
                if inst.features["class"]:
                    targets.append(1)
                else:
                    targets.append(0)
            #if imputer:
            #    imp = imputer
            #else:
            #    imp = Imputer(missing_values='NaN', strategy='median', axis=0)
            #    try:
            #        imp.fit(numpy.asarray(encoded_instances))
            #        cPickle.dump(imp, open('Imputer.dat', 'wb'), cPickle.HIGHEST_PROTOCOL)
            #    except:
            #        filename = str(uuid.uuid4())+'.dat'
            #        ofh = open(filename, "wb")
            #        cPickle.dump((encoded_instances, tags, copy.deepcopy(self.sets)), ofh, cPickle.HIGHEST_PROTOCOL)
            #        ofh.close()
            #        print "#########################################################################"
            #        print "###                                                                   ###"
            #        print "###                           ERROR OCCURRED                          ###"
            #        print "###  Instances File Dumped: %s  ###" % filename
            #        print "###                                                                   ###"
            #        print "#########################################################################"
            #        sys.exit(1)
            #encoding = (imp.transform(encoded_instances), numpy.asarray(targets, dtype=numpy.int), imp)
            encoding = (numpy.asarray(encoded_instances), numpy.asarray(targets, dtype=numpy.int))
        elif form == 'svm':
            encoding = ''
            for inst in instances:
                (y, x) = inst.get_instance(encoding=self, tags=tags, form='svm')
                encoding += str(y)+' '
                for key in x.keys():
                    encoding += str(key)+':'+str(x[key])+' '
                encoding += '\n'
        elif form == 'csv':
            header = 'class,'
            attrs = self.get_attrs(tags=tags, form='csv')
            for attr in attrs:
                header += str(attr)+','
            encoding = header[:-1]+'\n'
            for inst in instances:
                (y, x) = inst.get_instance(encoding=self, tags=tags, form='csv')
                encoding += str(y)+','+x
                encoding += '\n'
        return encoding

    def get_dataset(self, ratio=None, tags=None, imp=None):
        (desc, feature_names) = self.get_features(tags=tags, form='scikit')
        # (instances, targets, imp) = self.get_encoding(ratio=ratio, tags=tags, form='scikit', imputer=imp)
        (instances, targets) = self.get_encoding(ratio=ratio, tags=tags, form='scikit')
        target_names = ['0', '1']

        # return (Bunch(data=instances, feature_names=feature_names, DESC=desc, target=targets, target_names=target_names), imp)
        return residue.Bunch(data=instances, feature_names=feature_names, DESC=desc, target=targets, target_names=target_names)

    def write_arff(self, fh, ratio=None, tags=None):
        fh.write(self.get_features(tags=tags, form='arff'))
        fh.write(self.get_encoding(ratio=ratio, tags=tags, form='arff'))

    def write_csv(self, fh, ratio=None, tags=None):
        # fh.write(self.get_features(tags=tags, form='csv'))
        fh.write(self.get_encoding(ratio=ratio, tags=tags, form='csv'))
