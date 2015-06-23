__author__ = 'asyavuz'
#
# Warning! This file is a working draft. Complete version will be available till June 24, 2015.
#
import argparse
import cPickle
from Bio import SeqIO
from sequence import *

parser = argparse.ArgumentParser(description='Predicting neddylation target sites from protein sequences.')
parser.add_argument("FASTA", help="Input FASTA file")
parser.add_argument('pssm', help='Input PSSM file')
parser.add_argument('iupred', help='Input IUPred output file')
parser.add_argument('--thr', '-t', action='store', dest='thr', help='Prediction threshold (default: 0)')
parser.add_argument('--version', action='version', version='NeddyPreddy Standalone Version v0.1')
args = parser.parse_args()

if not args.thr:
    thr = 0.0
else:
    thr = float(args.thr)
(clf, scaler, metadata) = cPickle.load(open('NeddyPreddyModel.dat', 'rb'))

sequences = []
for seq_record in SeqIO.parse(args.FASTA, 'fasta'):

    seq = Seq(identifier=seq_record.id, sequence=str(seq_record.seq))
    seq.add_iupred(args.iupred)
    seq.add_pssm(args.pssm)
    seq.set_aa_freq()
    #print "Completed %s!" % (uniprot)
    sequences.append(seq)

residues = []
print "Identifier\tResidue\tClass\tDecision Value\tProbability"
for seq in sequences:
    for resid in seq:
        if resid == None:
            continue
        if resid.aa == 'K':
            sets = resid.set_selected(window_size=21)
            residues.append(resid)
            clf_ready = resid.get_instance(metadata['feats'], scaler).data

            decisionval = clf.decision_function(clf_ready)[0][0]

            dec_val = clf.decision_function(clf_ready)[0][0]
            prob = clf.predict_proba(clf_ready)[0][1]
            if dec_val>=thr:
                print "%s\tK%d\tNeddylated\t%.2f\t%.3f" % (seq.identifier, resid.no, dec_val, prob)
