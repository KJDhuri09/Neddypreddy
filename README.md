NeddyPreddy
================================
NeddyPreddy is a neddylation site predictor that uses SVM to predict neddylation sites from protein sequences. 

**[NEW]** NeddyPreddy web server is available at http://neddypreddy.sabanciuniv.edu (pending publication)

Installation
-------------
NeddyPreddy uses various python libraries, as well as BLAST+ executables and nr dataset. If you already have PSSM files of your proteins, you can skip BLAST+ installation.

You can get NeddyPreddy with 
```text
git clone https://asyavuz@bitbucket.org/asyavuz/neddypreddy.git
``` 
command.

Please refer individual websites of requirements for their installation instructions. 

### Requirements ###
* [BLAST+](http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [BioPython](http://biopython.org/)
* [scikit-learn](http://scikit-learn.org)
* [IUPred](http://iupred.enzim.hu/Downloads.php) *(for local disorder prediction)*
* [requests](http://docs.python-requests.org/en/master/) *(for online disorder prediction)*
* [BeautifulSoup4](https://www.crummy.com/software/BeautifulSoup/) *(for online disorder prediction)*

* * *

Usage
-----
The most basic use case would be as follows:

~~~~
python runneddypreddy.py my_protein.fasta
~~~~

For additional parameters please type:

~~~~
python runneddypreddy.py --help
~~~~

* * *

License
-------
NeddyPreddy uses the MIT License (see LICENCE.md). Please open an issue in this page if you have any questions.

* * *

Reference
---------
If you would like to use NeddyPreddy in your publications, please consider citing:
>  Yavuz AS, Sozer NB, Sezerman OS. "Prediction of neddylation sites from protein sequences 
>  and sequence-derived properties". BMC Bioinformatics, submitted.
