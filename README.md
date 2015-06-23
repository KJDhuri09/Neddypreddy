NeddyPreddy Standalone Version
------------------------------
You can clone this repository to run NeddyPreddy standalone version. Please beware that this is just a working draft at the moment. We will finalise the source code on June 24th. 

In order to run NeddyPreddy, you need to perform IUPred prediction and PSSM generation separately at the moment. 

Documentation will also be available on June 24, 2015.
Sorry for the inconvenience.

## Sample Usage ##

```
#!python

python runneddypreddy.py my_protein.fasta my_protein.pssm my_protein.iupred --thr 0.3
```

## Requirements ##
* BioPython: [http://biopython.org/](http://biopython.org/)
* scikit-learn: [http://scikit-learn.org](http://scikit-learn.org)