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
import logging
import subprocess, shlex
import tempfile
from multiprocessing import cpu_count
logging.basicConfig(filename='tools.log', level=logging.DEBUG)

MAX_CORE_PER_JOB = cpu_count()  # Change this value to 1 if you want psiblast to use single core only!

def obtain_pssm(sequence):
    '''
    Calculates PSSM with PSIBLAST. This function requires functional PSIBLAST installation and nr database.
    PSIBLAST should be reachable via PATH environment variable.

    :param sequence: Sequence object
    :return: Returns PSSM file name if PSSM calculation is successfully completed, False if calculation is failed.
    '''

    psiblast_exec = "psiblast"
    logger = logging.getLogger('psiblast')

    tmp = tempfile.NamedTemporaryFile(prefix='nptmp', suffix='.fasta', delete=False)
    seqfile = tmp.name
    tmp.write(sequence.get_fasta())
    tmp.close()
    logger.info("Running PSIBLAST for %s..." % sequence.identifier)

    command = psiblast_exec+" -num_threads %d -db nr -num_iterations 3 -inclusion_ethresh 1e-5 " % MAX_CORE_PER_JOB
    command += "-query "+seqfile+" -out_ascii_pssm "+seqfile+".pssm"
    logger.debug("Running command: %s" % command)
    args = shlex.split(command)
    try:
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except OSError:
        logger.error("Cannot find psiblast installation. Please make sure you installed BLAST+ and set PATH variable"
                     "correctly.")
        sys.stderr.write("Cannot find psiblast installation. Please make sure you installed BLAST+ and set "
                         "PATH variable correctly.\n")
        return False
    process.wait()
    stdout, stderr = process.communicate()

    if stdout:
        logger.info(stdout)
    if stderr:
        logger.error(stderr)

    logger.info("PSIBLAST run of "+seqfile+" finished with code "+str(process.returncode)+".")
    if int(process.returncode) == 0:
        os.remove(seqfile)
        return seqfile+".pssm"
    else:
        return False

def predict_disorder(sequence):
    '''
    Predicts disorder using IUPred. Tries to use local installation first, if fails, tries to use IUPred web server.
    Requires mechanize and BeautifulSoup modules.
    :param sequence: Sequence object
    :return: Returns IUPred prediction result file if prediction is successfully completed, False if prediction is failed.
    '''
    type = 'local'

    logger = logging.getLogger('iupred')

    tmp = tempfile.NamedTemporaryFile(prefix='nptmp', suffix='.fasta', delete=False)
    seqfile = tmp.name
    tmp.write(sequence.get_fasta())
    tmp.close()

    try:
        iupred = os.environ["IUPred_PATH"]+'/iupred'
    except KeyError:
        logger.info("IUPred environment variable is not set. Trying to access remote server...")
        sys.stderr.write("IUPred environment variable is not set. Trying to access remote server...\n")
        type = 'remote'

    out_file = seqfile+'.iupred'
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

        logger.info("IUPred run of "+sequence.identifier+" finished with code "+str(process.returncode)+".")
        if process.returncode != 0:
            logger.error("IUPred run was unsuccessful with return code %s." % str(process.returncode))
            return False
    else:
        try:
            import BeautifulSoup
        except ImportError:
            logger.error("Remote IUPred usage requires BeautifulSoup module, which is not installed in your system.")
            sys.stderr.write("Remote IUPred usage requires BeautifulSoup module, which is not installed in your system.\n")
            return False

        try:
            import mechanize
        except ImportError:
            logger.error("Remote IUPred usage requires mechanize module, which is not installed in your system.")
            sys.stderr.write("Remote IUPred usage requires mechanize module, which is not installed in your system.\n")
            return False

        logger.info("Loading IUPred server for %s." % sequence.identifier)
        response = mechanize.urlopen("http://iupred.enzim.hu/")
        forms = mechanize.ParseResponse(response, backwards_compat=False)
        try:
            form = forms[0]
            form['seq'] = sequence.sequence
            form['type'] = ['long']
            form['output'] = ['data']
            request2 = form.click()
            response2 = mechanize.urlopen(request2)
            content = response2.read()
        except:
            logger.error("Problem in IUPred remote server! Exception: %s" % (sys.exc_info()[0]))
            return False

        soup = BeautifulSoup.BeautifulSoup(content)
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
