
------------------------
Summary
------------------------
pplacer-XR (eXtra Range) is an extension for the maximum likelihood phylogenetic placement method pplacer, which allows pplacer to run on ultra-large reference trees. It is a python program that can be run on **Linux and MacOS**


------------------------
Requirements
------------------------
1. Python: Version >= 3.0
2. DendroPy: Version 4.4.0
3. pplacer: Version 1.1.alpha19


----------------------------------
Input & Output Specification
----------------------------------
The input parameters are as follows:

    model_name : DNA substitution model such as GTR or JC69
    statistics_file : statistics file produced by RAxML containing the substitution rates
    tree_path : refernce tree file path (in newick format) (parameter T)
    output_directory : directory to where output file pplacer-XR.jplace is written
    alnment : fasta format file containing the multiple sequence alignment (MSA) of reference and query sequences (parameter A)
    queries : txt file containing a list of query sequence labels (one per line) (parameter q)
    subtree_size : maximum size of the subtree for placement with pplacer (2000 is recommended) (parameter B) 
    random_number : arbitrary index of the tmp file created to contain the Newick subtree and its corresponding MSA

Usage:
Please run from the directory containing pplacer, available at https://github.com/matsen/pplacer.

    python3 pplacer-XR.py model_name statistic_file backbone_tree output_directory alignment queries subtree_size random_number
