NeddyPreddy
================================
NeddyPreddy is a neddylation site predictor that uses SVM to predict neddylation sites from protein sequences. 

Installation
-------------
NeddyPreddy uses various python libraries, as well as BLAST+ executables and nr database. If you already have PSSM files 
of your proteins, you can skip BLAST+ installation.

You can get NeddyPreddy with 
```text
git clone https://asyavuz@bitbucket.org/asyavuz/neddypreddy.git
``` 
command. 

NeddyPreddy is tested on Mac OSX (Yosemite) and Linux (CentOS 5.11, Ubuntu 14.04) systems. If you observe any problems, 
please report it at the [*Issues*](https://bitbucket.org/asyavuz/neddypreddy/issues) section.

Please refer individual websites of required software for their installation instructions. 

### Requirements ###
* [BLAST+](http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [BioPython](http://biopython.org/)
* [scikit-learn](http://scikit-learn.org) *(tested with 0.15.2 and 0.16.1)*
* [IUPred](http://iupred.enzim.hu/Downloads.php) *(for local disorder prediction)*
* [mechanize](http://wwwsearch.sourceforge.net/mechanize/) *(for online disorder prediction)*
* [BeautifulSoup](http://www.crummy.com/software/BeautifulSoup/) *(for online disorder prediction)*

* * *

Usage
-----
The most basic use case would be as following:

~~~~
python runneddypreddy.py my_protein.fasta
~~~~

For additional parameters please type:

~~~~
python runneddypreddy.py --help
~~~~

* * *

Licence
-------
NeddyPreddy uses the MIT License (see [LICENCE.md](https://bitbucket.org/asyavuz/neddypreddy/raw/master/LICENSE.md)). Please open an issue in this page if you have any questions.

* * *

Reference
---------
If you would like to use NeddyPreddy in your publications, please consider citing:
>  Yavuz AS, Sozer NB, Sezerman OS. "Prediction of neddylation sites from protein sequences 
>  and sequence-derived properties". *BMC Bioinformatics*, submitted.
