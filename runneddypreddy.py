__author__ = 'asyavuz'
#
# runneddypreddy.py
# Predicts neddylation sites from protein sequences.
#
# Author: Ahmet Sinan Yavuz <asinanyavuz@sabanciuniv.edu>
#
# Please obtain latest version and documentation from https://bitbucket.org/asyavuz/neddypreddy
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
import os
import sys
import argparse
import cPickle
from Bio import SeqIO
from sequence import *
from tools import *

# Command line argument parser
parser = argparse.ArgumentParser(description='Predicting neddylation target sites from protein sequences.')
parser.add_argument("FASTA", help="Input protein sequence file (FASTA-formatted)")
parser.add_argument('--fmt', '-f', action='store', dest='fmt', type=str, default='tsv',
                    help='Output format (available: CSV or TSV, default: TSV)')
parser.add_argument('--out', '-o', action='store', dest='out', help='Output filename (default: standard output)')
parser.add_argument('--pssm', '-p', action='store', dest='pssm', help='Input PSSM file (only works for single protein '
                                                                      'prediction)')
parser.add_argument('--dis', '-d', action='store', dest='dis', help='Input IUPred output file (only works for single '
                                                                    'protein prediction)')
parser.add_argument('--thr', '-t', action='store', dest='thr',
                    help='Prediction threshold (accepts: floating number between -2.0 and 2.0, or \'low\',''medium\', '
                         '\'high\', default: 0.0)')
parser.add_argument('--version', '-v', action='version', version='NeddyPreddy Standalone Version v0.2')
args = parser.parse_args()

if not args.thr:  # Set default threshold, if its not specified by the user
    thr = (0.0, 'medium')
else:
    if args.thr.lower() == 'low':  # Convert user provided keywords to actual thresholds
        thr = (-0.9, 'low')
    elif args.thr.lower() == 'medium':
        thr = (0.0, 'medium')
    elif args.thr.lower() == 'high':
        thr = (1.0, 'high')
    else:
        thr = (float(args.thr), 'custom')

# File check
if not os.path.exists(args.FASTA):
    sys.stderr.write("Cannot find input protein sequence file: %s\n" % args.FASTA)
    sys.exit(1)
if args.pssm and not os.path.exists(args.pssm):
    sys.stderr.write("Cannot find PSSM file: %s\n" % args.pssm)
    sys.exit(1)
if args.dis and not os.path.exists(args.dis):
    sys.stderr.write("Cannot find disorder file: %s\n" % args.dis)
    sys.exit(1)

print "NeddyPreddy - neddylation site prediction"
print "Loading input file: %s" % args.FASTA
# Parse sequences
sequence_records = list(SeqIO.parse(args.FASTA, 'fasta'))

if len(sequence_records) > 1 and (args.pssm or args.dis):
    sys.stderr.write("Cannot use PSSM and disorder predictions in multiple sequence prediction! Ignoring files...")
    args.dis, args.pssm = '', ''

sequences = []
for seq_record in sequence_records:
    seq = Seq(identifier=seq_record.id, sequence=str(seq_record.seq))

    if args.dis:
        seq.add_iupred(args.dis)
    else:
        print "Performing IUPred prediction for %s..." % seq.identifier
        tmpfile = predict_disorder(seq)  # Predict disorder from scratch (requires IUPred or internet connection)
        if tmpfile:
            seq.add_iupred(tmpfile)
        else:
            sys.stderr.write("Error in disorder prediction. Please check tools.log.\n")
            sys.exit(1)

    if args.pssm:
        seq.add_pssm(args.pssm)
    else:
        print "Calculating PSSM for %s..." % seq.identifier
        tmpfile = obtain_pssm(seq)  # Obtain PSSM from protein sequence (requires PSIBLAST and nr database)
        if tmpfile:
            seq.add_pssm(tmpfile)
        else:
            sys.stderr.write("Error in PSSM calculation. Please check tools.log.\n")
            sys.exit(1)

    seq.set_aa_freq()
    sequences.append(seq)

# Generate classifier object
(clf, scaler, metadata) = get_classifier()
if not clf:
    sys.stderr.write("An issue occurred during classifier training. Please check tools.log.\n")
    sys.exit(1)

if os.path.exists('tools.log'):  # If tools have been used, removing log here, as you wouldn't need it after this point
    os.remove('tools.log')

prediction_results =[]
for seq in sequences:
        for resid in seq:
            if resid == None:  # Skip placeholder aminoacid
                continue
            if resid.aa == 'K':
                resid.set_selected()

                clf_ready = resid.get_instance(metadata['feats'], scaler).data

                decfunout = clf.decision_function(clf_ready)
                if decfunout[0] is list: # scikit-learn returns decision value in a nested list in 0.15, in a list in 0.16
                    dec_val = decfunout[0][0]
                else:
                    dec_val = decfunout[0]

                prob = clf.predict_proba(clf_ready)[0][1]

                if dec_val >= thr[0]:
                    prediction_results.append((seq.identifier, 'K'+str(resid.no), 'Neddylated',
                                               '%.2f' % dec_val, '%.2f' % prob))

if args.out:  # Set output stream
    ostream = open(args.out, 'w')
else:
    ostream = sys.__stdout__

if args.fmt.lower() == 'csv':  # Set separator for results
    sep = ','
else:
    sep = '\t'

# Output results
ostream.write("# NeddyPreddy - Prediction Results\n")
ostream.write("# Model Version: %s\n" % str(metadata['Date']))
ostream.write("# Threshold: %.2f (%s)\n" % thr)
ostream.write(sep.join(["Identifier", "Residue", "Class", "Decision Value", "Probability"])+'\n')
for res in prediction_results:
    ostream.write(sep.join(res)+'\n')

if args.out:
    print "Completed."
