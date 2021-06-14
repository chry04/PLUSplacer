

------------------------
Summary
------------------------
pplacer-XR (eXtra Range) is an extension for the maximum likelihood phylogenetic placement method pplacer, which allows pplacer to run on ultra-large reference trees. 

Similarly EPA-ng-XR (eXtra Range) is an extention of EPA-ng, allowing it to run on ultra-large reference trees.

They are both python programs that can be run on **Linux and MacOS**


------------------------
Requirements
------------------------
1. Python: Version >= 3.0
2. Treeswift
3. pplacer: Version 1.1.alpha19 (for pplacer-XR)
   or 
   EPA-ng: Version 0.3.8 (for EPA-ng-XR)


----------------------------------
Input & Output Specification
----------------------------------
The input parameters for pplacer-XR and EPA-ng-XR are as follows:
    
    Required arguments: 
    -i, --info : statistics file produced by RAxML version 7 for pplacer-XR or RAxML-ng for EPA-ng-XR containing the substitution rates
    -t, --tree : refernce tree file path (in newick format) (parameter T)
    -d, --outdir : directory to where output file pplacer-XR.jplace is written
    -a, --alignment : fasta format file containing the multiple sequence alignment (MSA) of reference and query sequences (parameter A) (query sequences optional if aligned queries are in a separate file)

    Optional arguments:
    -m, --model : DNA substitution model such as GTR for nucleotides
    -q, --qalignment : file containing a list of aligned query sequences in fasta format (needed when not contained in file with reference MSA)
    -b, --subtreesize : maximum size of the subtree for placement (with pplacer 2000 is recommended, and with EPA-ng 10,000 is recommended) (parameter B) 
    -s, --subtreetype : options for collecting nodes for the subtree - d (default) for edge weighted distances, n for node distances, h for hamming distances
    -n, --tmpfilenbr : number for working directories in the current file path (to allow to run multiple instances concurently)
    -f, --fragmentflag : boolean, True is you wish to mask leading and trailing gaps in each query sequence when finding closest sister taxon.
    

Usage:
Please run from a directory containing the respective phylogenetic placement method. This would be pplacer (available at https://github.com/matsen/pplacer) for pplacer-XR or EPA-ng (available at https://github.com/Pbdas/epa-ng) for EPA-ng-XR.

python3 pplacer-XR.py -i INFO -t TREE -d OUTDIR -a ALIGNMENT

python3 EPA-ng-XR.py -i INFO -t TREE -d OUTDIR -a ALIGNMENT    
