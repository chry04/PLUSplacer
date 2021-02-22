"""
This file is contains code to be used alongside pplacer as described in the upcoming ALCOB conference paper :

Wedell, E., Cai, Y., Warnow, T. (2021). Scalable and Accurate Phylogenetic Placementusing pplacer-XR.

Copyright (c) 2021 pplacer-XR Developers
Yirong Cai <yirongc2@illinois.edu>
Eleanor Wedell <ewedell2@illinois.edu>
All rights reserved.

Licence: MIT Licence, 
see https://opensource.org/licenses/MIT

pplacer can be found at: 
https://github.com/matsen/pplacer

***must be run from the directory containing pplacer***
"""

import sys
import os
import utils
import dendropy
import shutil
import json
from dendropy.calculate import treecompare

if __name__ == "__main__":
    model = str(sys.argv[1])
    info = str(sys.argv[2])
    tree_path = str(sys.argv[3])
    output = str(sys.argv[4]) 
    aln = str(sys.argv[5]) 
    query = str(sys.argv[6]) 
    n = int(sys.argv[7]) 
    run = int(sys.argv[8]) 
    
    """
    Parameters
    ----------
    model : DNA substitution model such as GTR or JC69
    info : statistics file produced by RAxML containing the substitution rates
    tree_path : refernce tree file path (in newick format) (parameter T)
    output : directory to where output file pplacer-XR.jplace is written
    aln : fasta format file containing the multiple sequence alignment (MSA) of reference and query sequences (parameter A)
    query : txt file containing a list of query sequence labels (one per line) (parameter q)
    n : maximum size of the subtree for placement with pplacer (parameter B)
    run : arbitrary index of the tmp file created to contain the Newick subtree and its corresponding MSA
    
    Output
    ------
    pplacer-XR.jplace : placement file in jplace format
    """

    # read msa and reference tree
    aln_dict = utils.read_data(aln)
    ref_dict, q_dict = utils.seperate(aln_dict, query)
    tree = dendropy.Tree.get(path=tree_path, schema="newick", rooting='force-unrooted')
    
    tree.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    namespace = tree.taxon_namespace
    jplace = dict()
    placements = []
    
    # create edge tokens to be used with jplace output
    counter = 0
    for e in tree.postorder_edge_iter():
        e.label = counter
        counter += 1
    jplace["tree"] = tree.as_string(schema="newick", suppress_edge_lengths=False, \
        edge_label_compose_fn=utils.edge_labeling, suppress_rooting=True)[:-1]

    # create output directory and tmp file
    files = []
    try: 
        os.mkdir("tmp{}".format(run))
    except:
    	print("tmp directory already exists")    
    try:
        os.mkdir(output)
    except:
    	print("path directory already exists")
        
    # place each query sequence
    for name, seq in q_dict.items():
        tmp_tree = "../scratch/tmp{}/tree_".format(run)+name
        tmp_aln = "../scratch/tmp{}/aln".format(run)+name+".fa"
        tmp_output = output+"/"+name+".jplace"
        
        # skip if the jplace file placing that sequence already exists
        if not os.path.exists(tmp_output):
            #finds closest sister taxon l
            y = utils.find_y(seq, ref_dict)
            
            #create subtree containing closest sister taxon l
            nodes = utils.subtree_nodes(tree, y, n)
            subtree = tree.extract_tree_with_taxa(nodes)
            subtree.collapse_basal_bifurcation(set_as_unrooted_tree=True)
            subtree.write(path=tmp_tree, schema="newick", suppress_rooting=True)

            # write subtree MSA and aligned query sequence to tmp file
            f = open(tmp_aln, "w")
            f.write(">"+name)
            f.write("\n")
            f.write(seq+"\n")
            for ns in nodes:
                f.write(">"+ns.label+"\n")
                f.write(ref_dict[ns.label])
                f.write("\n")

            f.close()

            # run pplacer from directory containing pplacer binaries
            os.system("./pplacer -m {} -s {} -t {} --keep-at-most 1 -o {} {}".format(model, info, tmp_tree, tmp_output, tmp_aln))

        # load the jplace file and find placements in the original backbone tree
        place_file = open(tmp_output, 'r')
        place_json = json.load(place_file)

        if len(place_json["placements"]) > 0:
            
            added_tree = dendropy.Tree.get(data=place_json["tree"], 
            				    schema="newick", 
                                            is_parse_jplace_tokens = True,
        				    taxon_namespace=namespace)
        
            tmp_place = place_json["placements"][0]
            
            for i in range(len(tmp_place["p"])):
                edge_num = tmp_place["p"][i][1] # edge number in subtree
                edge_distal = tmp_place["p"][i][0] # distal length from parent node
                
                # find placement edge according to edge number
                place_edge = added_tree.edges(\
                        filter_fn=(lambda x: str(x.edge_number)==str(edge_num)))[0]
                right_n = place_edge.head_node
                left_n = right_n.parent_node
                
                # obtain a path from leaf left to leaf right containing placement edge through the subtree
                left, path_l = utils.find_closest(left_n, {left_n, right_n})
                right, path_r = utils.find_closest(right_n, {left_n, right_n})
                
                # obtain the corresponding path in backbone tree
                left = tree.find_node_with_taxon_label(left.taxon.label)
                right = tree.find_node_with_taxon_label(right.taxon.label)
                _, path = utils.find_closest(left, {left}, y=right)

                # find the length of placement along the path from leaf left to leaf right in subtree
                length = sum([x.length for x in path_l])+edge_distal 
                
                # find the target placement edge in backbone tree 
                target_edge = path[-1]
                for j in range(len(path)):
                    length -= path[j].length
                    if length < 0:
                        target_edge = path[j]
                        break
                                               
                tmp_place["p"][i][0] = 0

                tmp_place["p"][i][0] = target_edge.length+length
                tmp_place["p"][i][1] = target_edge.label
                                
            # append the placement to the output jplace
            placements.append(tmp_place.copy())

        place_file.close()

    # build jplace file
    jplace["placements"] = placements
    jplace["metadata"] = {"invocation": " ".join(sys.argv)}
    jplace["version"] = 3
    jplace["fields"] = ["distal_length", "edge_num", "like_weight_ratio", \
            "likelihood", "pendant_length"]

    output = open(output+"/pplacer-XR.jplace", 'w')
    json.dump(jplace, output, sort_keys=True , indent=4)
    output.close()
    
    #remove tmp directory 
    shutil.rmtree("tmp{}".format(run))
