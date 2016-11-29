__author__ = 'asyavuz'
#
# tools.py
# Helper functions for NeddyPreddy Standalone Version
#
# Author: Ahmet Sinan Yavuz <asinanyavuz@sabanciuniv.edu>
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
import json
import logging
import subprocess, shlex
import tempfile
import numpy as np
import sklearn as sk
from multiprocessing import cpu_count
import sklearn
from sklearn import svm, preprocessing

logging.basicConfig(filename='tools.log', level=logging.DEBUG)

MAX_CORE_PER_JOB = cpu_count()  # Change this value to 1 if you want psiblast to use single core only!


class Bunch(dict):
    """Container object for datasets: dictionary-like object that
       exposes its keys as attributes."""

    def __init__(self, **kwargs):
        dict.__init__(self, kwargs)
        self.__dict__ = self


def obtain_pssm(sequence):
    '''
    Calculates PSSM with PSIBLAST. This function requires functional PSIBLAST installation and nr database.
    PSIBLAST should be reachable via PATH environment variable.

    :param sequence: Sequence object
    :return: Returns PSSM file name if PSSM calculation is successfully completed, False if calculation is failed.
    '''

    psiblast_exec = "psiblast"
    logger = logging.getLogger('psiblast')

    tmp = tempfile.NamedTemporaryFile(mode="w", prefix='nptmp', suffix='.fasta', delete=False)
    seqfile = tmp.name
    tmp.write(sequence.get_fasta())
    tmp.close()

    tmp_stdout = tempfile.TemporaryFile()
    tmp_stderr = tempfile.TemporaryFile()
    logger.info("Running PSIBLAST for %s..." % sequence.identifier)

    command = psiblast_exec + " -num_threads %d -db nr -num_iterations 3 -inclusion_ethresh 1e-5 " % MAX_CORE_PER_JOB
    command += "-query " + seqfile + " -out_ascii_pssm " + seqfile + ".pssm"
    logger.debug("Running command: %s" % command)
    args = shlex.split(command)
    try:
        process = subprocess.Popen(args, stdout=tmp_stdout, stderr=tmp_stderr)
    except OSError:
        logger.error("Cannot find psiblast installation. Please make sure you installed BLAST+ and set PATH variable"
                     "correctly.")
        sys.stderr.write("Cannot find psiblast installation. Please make sure you installed BLAST+ and set "
                         "PATH variable correctly.\n")
        return False
    process.wait()
    # stdout, stderr = process.communicate()

    # if stdout:
    #    logger.info(stdout)
    # if stderr:
    #    logger.error(stderr)

    logger.info("PSIBLAST run of " + seqfile + " finished with code " + str(process.returncode) + ".")
    if int(process.returncode) == 0:
        os.remove(seqfile)
        return seqfile + ".pssm"
    else:
        return False


def predict_disorder(sequence):
    '''
    Predicts disorder using IUPred. Tries to use local installation first, if fails, tries to use IUPred web server.
    Requires mechanize and BeautifulSoup4 modules.
    :param sequence: Sequence object
    :return: Returns IUPred prediction result file if prediction is successfully completed, False if prediction is failed.
    '''
    type = 'local'

    logger = logging.getLogger('iupred')

    tmp = tempfile.NamedTemporaryFile(mode="w", prefix='nptmp', suffix='.fasta', delete=False)
    seqfile = tmp.name
    tmp.write(sequence.get_fasta())
    tmp.close()

    try:
        iupred = os.environ["IUPred_PATH"] + '/iupred'
    except KeyError:
        logger.info("IUPred environment variable is not set. Trying to access remote server...")
        sys.stderr.write("IUPred environment variable is not set. Trying to access remote server...\n")
        type = 'remote'

    out_file = seqfile + '.iupred'
    if not type == 'remote':
        ofh = open(out_file, 'w')
        command = "%s %s long" % (iupred, seqfile)
        logger.info("Running IUPred with command: %s" % command)
        logger.info("Expected output file: %s" % out_file)
        args = shlex.split(command)
        process = subprocess.Popen(args, stdout=ofh, stderr=subprocess.PIPE)
        process.wait()
        ofh.close()

        stdout, stderr = process.communicate()

        if stdout:
            logger.info(stdout)
        if stderr:
            logger.error(stderr)

        logger.info("IUPred run of " + sequence.identifier + " finished with code " + str(process.returncode) + ".")
        if process.returncode != 0:
            logger.error("IUPred run was unsuccessful with return code %s." % str(process.returncode))
            return False
    else:
        try:
            from bs4 import BeautifulSoup
        except ImportError:
            logger.error("Remote IUPred usage requires BeautifulSoup4 module, which is not installed in your system.")
            sys.stderr.write(
                "Remote IUPred usage requires BeautifulSoup4 module, which is not installed in your system.\n")
            return False

        try:
            import requests
        except ImportError:
            logger.error("Remote IUPred usage requires requests module, which is not installed in your system.")
            sys.stderr.write("Remote IUPred usage requires requests module, which is not installed in your system.\n")
            return False

        logger.info("Loading IUPred server for %s." % sequence.identifier)
        try:
            payload = {"seq": sequence.sequence, "type": "long", "output": "data", "sp": ""}
            response = requests.post("http://iupred.enzim.hu/pred.php", data=payload)
            content = response.text
        except Exception as e:
            logger.error("Problem in IUPred remote server! Exception: %s" % str(e))
            return False

        soup = BeautifulSoup(content, "html.parser")
        table = soup.find('table', attrs={'width': 420})

        ofh = open(out_file, 'w')

        rows = table.findAll('tr')
        for i, tr in enumerate(rows):
            if i == 1:
                continue

            row_content = []
            cols = tr.findAll('td')
            for j, td in enumerate(cols):
                text = ''.join(td.find(text=True))
                if (i == 0 and j == 0) and str(text) != "Disorder prediction score":
                    break
                elif i >= 2:
                    row_content.append(str(text.strip()))
            if row_content:
                ofh.write("%s\t%s\t%s\n" % (row_content[0], row_content[1], row_content[2]))
        ofh.close()
        logger.info("Completed prediction for %s, file saved to: %s." % (sequence.identifier, out_file))

    os.remove(seqfile)
    return out_file


def get_classifier(dname):
    '''
    This function fits a classifier model to the training data. Classifier is trained in the each run to ensure
    compability between different scikit-learn versions. Serialized classifier objects fail to function in different
    versions of scikit-learn.
    :return: classifier and scaler object, as well as metadata of classifier.
    '''
    logger = logging.getLogger('classifier')

    if sklearn.__version__ >= "0.18":
        from sklearn.model_selection import cross_val_score
    else:
        from sklearn.cross_validation import cross_val_score

    train_dataset = Bunch()
    train_dataset.data = np.load(dname + "/data/TS.npy")
    train_dataset.target = np.load(dname + "/data/TT.npy")
    with open(dname + "/data/Metadata.json") as fh:
        metadata = json.load(fh)

    try:
        scaler = preprocessing.MinMaxScaler().fit(train_dataset.data)
    except Exception as e:
        logger.error("Problem in scaler fitting. Exception: %s" % str(e))
        return False, False, False

    try:
        # Fix class weights, json does not allow numeric keys, scikit does not like string keys
        metadata["class_weight"][0] = metadata["class_weight"]["0"]
        metadata["class_weight"][1] = metadata["class_weight"]["1"]
        del metadata["class_weight"]["0"]
        del metadata["class_weight"]["1"]

        clf = svm.SVC(C=metadata['C'], kernel='rbf', gamma=metadata['gamma'], class_weight=metadata['class_weight'],
                      probability=True)

        # Check if cross-validation score is consistent with published method (for compatibility)
        scores = cross_val_score(clf, scaler.transform(train_dataset.data), train_dataset.target, cv=5)
        assert scores.mean() - (scores.std() * 2) <= 0.91 <= scores.mean() + (scores.std() * 2)

        clf.fit(scaler.transform(train_dataset.data), train_dataset.target)
    except Exception as e:
        logger.error("Problem in training the classifier. Exception: %s" % str(e))
        sys.stderr.write("Problem in training the NeddyPreddy classifier.\n")
        return False, False, False

    return clf, scaler, metadata
