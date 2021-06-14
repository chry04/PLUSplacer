"""
This file is contains code to be used alongside pplacer as described in the upcoming ALCOB conference paper :

Wedell, E., Cai, Y., Warnow, T. (2021). Scalable and Accurate Phylogenetic Placementusing pplacer-XR.

Copyright (c) 2021 EPA-ng-XR Developers
Yirong Cai <yirongc2@illinois.edu>
Eleanor Wedell <ewedell2@illinois.edu>
All rights reserved.

Licence: MIT Licence, 
see https://opensource.org/licenses/MIT

EPA-ng can be found at: 
https://github.com/Pbdas/epa-ng

***must be run from the directory containing EPA-ng***
"""


import sys
import os
import utils
import shutil
import json
import time
import argparse
import treeswift

def main(args):
    tree_path = args.tree
    output = args.outdir
    outFile = args.output
    aln = args.alignment
    n = args.subtreesize
    run = args.tmpfilenbr
    subtree_flag = args.subtreetype
    frag_flag = args.fragmentflag
    q_aln = args.qalignment
    model = args.model
    info = args.info
    
    # read msa and reference tree
    t0 = time.perf_counter()
    tree = treeswift.read_tree_newick(tree_path)
    tree.deroot()
    leaf_dict = tree.label_to_node(selection='leaves')

    if q_aln != "":
        ref_dict = utils.read_data(aln)
        q_dict = utils.read_data(q_aln)
    else:
        aln_dict = utils.read_data(aln)
        ref_dict, q_dict = utils.seperate(aln_dict, leaf_dict)

    jplace = dict()
    placements = []

    # create edge tokens to be used with jplace output
    utils.add_edge_nbrs(tree)
    jplace["tree"] = utils.newick_edge_tokens(tree)
    print ('{} seconds loading data'.format(time.perf_counter() - t0))

    files = []
    try:
        os.mkdir("tmp{}".format(run))
    except OSError as error:
    	pass
    try:
        os.mkdir(output)
    except OSError as error:
    	pass


    # place each query sequence
    for name, seq in q_dict.items():

        tmp_tree = "tmp{}/tree_".format(run) + name
        tmp_aln = "tmp{}/aln".format(run) + name + ".fa"
        tmp_qaln = "tmp{}/qaln".format(run) + name + "q.fa"
        tmp_output = "tmp{}/{}".format(run, name) + "/epa_result.jplace"
        tmp_dir = "tmp{}/{}".format(run, name)
        try:
            os.mkdir(tmp_dir)
        except OSError as error:
            pass


        # finds closest sister taxon and subtree leaves
        if subtree_flag == 'h':
            y = utils.find_closest_hamming(seq, ref_dict, n, frag_flag)
            labels = []
            for taxon in y:
               labels.append(leaf_dict[taxon].get_label())
            print('Closest sister taxon found: {}'.format(y[0]))
        else:
            y = utils.find_closest_hamming(seq, ref_dict, 1, frag_flag)
            print ('Closest sister taxon found: {}'.format(y[0]))
            print ('{} seconds finding closest leaf'.format(time.perf_counter() - t0))
            node_y = leaf_dict[y[0]]
            if subtree_flag == 'n':
                labels = utils.subtree_nodes(tree, node_y, n)
            else:
                labels = utils.subtree_nodes_with_edge_length(tree, node_y, n)

        # write subtree MSA and aligned query sequence to tmp file
        f = open(tmp_aln, "w")
        fq = open(tmp_qaln, "w")
        fq.write(">"+name)
        fq.write("\n")
        fq.write(seq+"\n")
        for label in labels:
            label_list = label.split('%%',1)
            label = label_list[0]
            f.write(">"+label+"\n")
            f.write(ref_dict[label])
            f.write("\n")

        f.close()
        fq.close()
        

        subtree = tree.extract_tree_with(labels)

        subtree.deroot()
        utils.remove_edge_nbrs(subtree)

        subtree.write_tree_newick(tmp_tree, hide_rooted_prefix=True)

        print ('{} seconds extracting subtree'.format(time.perf_counter() - t0))
        # run EPA-ng-XR from directory containing EPA-ng binaries
        os.system("./epa-ng -m {} -t {} -w {} -s {} -q {} --redo -T 16".format(info, tmp_tree, tmp_dir, tmp_aln, tmp_qaln))

        print ('{} seconds running epa-ng'.format(time.perf_counter() - t0))

        # load the jplace file and find placements in the original backbone tree
        place_file = open(tmp_output, 'r')
        place_json = json.load(place_file)

        if len(place_json["placements"]) > 0:

            added_tree, edge_dict = utils.read_tree_newick_edge_tokens(place_json["tree"])

            tmp_place = place_json["placements"][0]
            for i in range(len(tmp_place["p"])):
                edge_num = tmp_place["p"][i][0] # edge number in subtree
                edge_distal = tmp_place["p"][i][3] # distal length from parent node

                # find placement edge according to edge number
                right_n = edge_dict[str(edge_num)]
                left_n = right_n.get_parent()

                # obtain a path from leaf left to leaf right containing placement edge through the subtree
                left, path_l = utils.find_closest(left_n, {left_n, right_n})
                right, path_r = utils.find_closest(right_n, {left_n, right_n})

                # obtain the corresponding path in backbone tree
                left = leaf_dict[left.get_label()]
                right = leaf_dict[right.get_label()]
                _, path = utils.find_closest(left, {left}, y=right)

                # find the length of placement along the path from leaf left to leaf right in subtree
                length = sum([x.get_edge_length() for x in path_l])+edge_distal

                # find the target placement edge in backbone tree 
                target_edge = path[-1]
                for j in range(len(path)):
                    length -= path[j].get_edge_length()
                    if length < 0:
                        target_edge = path[j]
                        break

                tmp_place["p"][i][0] = 0

                label = target_edge.get_label()
                [taxon, target_edge_nbr] = label.split('%%',1)
                tmp_place["p"][i][0] = target_edge.get_edge_length()+length
                tmp_place["p"][i][1] = int(target_edge_nbr)

            # append the placement to the output jplace
            placements.append(tmp_place.copy())

        place_file.close()

    # build jplace file
    jplace["placements"] = placements
    jplace["metadata"] = {"invocation": " ".join(sys.argv)}
    jplace["version"] = 3
    jplace["fields"] = ["distal_length", "edge_num", "like_weight_ratio", \
            "likelihood", "pendant_length"]


    output = open('{}/{}.jplace'.format(output,outFile), 'w')
    json.dump(jplace, output, sort_keys=True , indent=4)
    output.close()
    print ('{} seconds building jplace'.format(time.perf_counter() - t0))
    shutil.rmtree("tmp{}".format(run))
    
    
def parseArgs():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--info", type=str,
                        help="Path to model parameters", required=True, default=None)
    
    parser.add_argument("-t", "--tree", type=str,
                        help="Path to reference tree with estimated branch lengths", required=True, default=None)
    
    parser.add_argument("-d", "--outdir", type=str,
                        help="Directory path for output", required=True, default=None)
    
    parser.add_argument("-a", "--alignment", type=str,
                        help="Path for query and reference sequence alignment in fasta format", required=True, default=None)

    parser.add_argument("-o", "--output", type=str,
                        help="Output file name", required=False, default="pplacer-XR")
    
    parser.add_argument("-m", "--model", type=str,
                        help="Model used for edge distances",
                        required=False, default="GTR")

    parser.add_argument("-b", "--subtreesize", type=int,
                        help="Integer size of the subtree",
                        required=False, default=10000)
    
    parser.add_argument("-s", "--subtreetype", type=str,
                        help="d (default) for edge weighted distances, n for node distances, h for hamming distances",
                        required=False, default="d")
    
    parser.add_argument("-n","--tmpfilenbr", type=int,
                        help="tmp file number",
                        required=False, default=0)
    
    parser.add_argument("-q", "--qalignment", type=str,
                        help="Path to query sequence alignment in fasta format (ref alignment separate)",
                        required=False, default="")
    
    parser.add_argument("-f", "--fragmentflag", type=bool,
                        help="boolean, True if queries contain fragments",
                        required=False, default=False)

    parser.add_argument("-v", "--version", action="version", version="2.0.0", help="show the version number and exit")
                       
    return parser.parse_args()

if __name__ == "__main__":
    main(parseArgs())
