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

import concurrent.futures
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor

# my own preset - Chengze Shen
# default settings for tqdm progress bar style
tqdm_styles = {
        'desc': '\tRunning...', 'ascii': False,
        'ncols': 80, 
        #'disable': True,
        'colour': 'green',
        'mininterval': 0.5
        }

__root = os.path.dirname(os.path.realpath(__file__))
# find existing epa-ng install
epa_ng_bin = shutil.which('epa-ng')

global tree, leaf_dict, ref_dict, query_ham_dict 

def initializer(_tree, _leaf_dict, _ref_dict, _query_ham_dict):
    global tree, leaf_dict, ref_dict, query_ham_dict
    tree = _tree
    leaf_dict = _leaf_dict
    query_ham_dict = _query_ham_dict
    ref_dict = _ref_dict

def runner(args):
    t0 = time.perf_counter()
    output, run, qnum, name, seq, n, info, subtree_flag, num_threads = args
    tmp_dir = "{}/tmp{}/{}".format(output, run, qnum)
    tmp_tree = "{}/tmp{}/{}/tree".format(output, run, qnum)
    tmp_aln = "{}/tmp{}/{}/aln.fa".format(output, run, qnum)
    tmp_qaln = "{}/tmp{}/{}/qaln.fa".format(output, run, qnum)
    tmp_output = "{}/tmp{}/{}/epa_result.jplace".format(output, run, qnum)
    try:
        os.mkdir(tmp_dir)
    except OSError as error:
        return None, -1, -1

    # finds closest sister taxon and subtree leaves
    if subtree_flag == 'h':
        y = query_ham_dict[name]
        labels = []
        for taxon in y:
           labels.append(leaf_dict[taxon].get_label())
    else:
        y = query_ham_dict[name]
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

    utils.remove_edge_nbrs(subtree)

    subtree.write_tree_newick(tmp_tree, hide_rooted_prefix=True)

    t1 = time.perf_counter()
    time_extract_subtree = t1 - t0
    #print ('{} seconds extracting subtree'.format(time.perf_counter() - t0))
    # run EPA-ng-XR from directory containing EPA-ng binaries
    os.system("{} -m {} -t {} -w {} -s {} -q {} --redo -T {}".format(
        epa_ng_bin, info, tmp_tree, tmp_dir, tmp_aln, tmp_qaln, 1))

    t2 = time.perf_counter()
    time_epa_ng = t2 - t1
    #print ('{} seconds running epa-ng'.format(time.perf_counter() - t0))

    return tmp_output, time_extract_subtree, time_epa_ng

def main(args):
    tree_path = args.tree
    output = args.outdir
    outFile = args.output
    aln = args.alignment
    n = args.subtreesize
    run = args.tmpfilenbr
    subtree_flag = args.subtreetype
    fragment_flag = args.fragmentflag
    q_aln = args.qalignment
    model = args.model
    info = args.info
    num_threads = args.threads
    
    if epa_ng_bin == '':
        raise ValueError('epa-ng binary not found in PATH!')
    
    # read msa and reference tree
    t0 = time.perf_counter()
    tree = treeswift.read_tree_newick(tree_path)
    tree.deroot()
    leaf_dict = tree.label_to_node(selection='leaves')

    # clean the leaf keys so that ' or " are not present
    ori_keys = list(leaf_dict.keys())
    for key in ori_keys:
        _node = leaf_dict[key]
        new_key = key.replace('\'', '')
        new_key = new_key.replace('\"', '')
        leaf_dict.pop(key)
        leaf_dict[new_key] = _node

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

    try:
        os.mkdir("{}/tmp{}".format(output, run))
    except OSError as error:
    	pass
    try:
        os.mkdir(output)
    except OSError as error:
    	pass

    # find the nearest taxon for each query
    tmp_output = "{}/tmp{}/".format(output, run) + "closest.txt"
    
    if q_aln == "":
        q_aln = "{}/tmp{}/".format(output, run) + "qaln.fa"
        f = open(q_aln, "w")
        for label, seq in q_dict.items():
            f.write(">"+label+"\n")
            f.write(seq+"\n")
        f.close()
        
        aln = "{}/tmp{}/".format(output, run) + "aln.fa"
        f = open(aln, "w")
        for label, seq in ref_dict.items():
            f.write(">"+label+"\n")
            f.write(seq+"\n")
        f.close()
    
    if subtree_flag == 'h':
        nbr_closest = n
    else:
        nbr_closest = 1

    if fragment_flag:
        os.system("{}/fragment_hamming {} {} {} {} {} {}".format(
            __root, aln, len(ref_dict), q_aln, len(q_dict), tmp_output, nbr_closest))
    else:
        os.system("{}/hamming {} {} {} {} {} {}".format(
            __root, aln, len(ref_dict), q_aln, len(q_dict), tmp_output, nbr_closest))
        
    print ('{} seconds finding closest leaves'.format(time.perf_counter() - t0))

    query_ham_dict = dict()
    f = open(tmp_output)
    for line in f: 
        line = line.strip()
        y = line.split(',') 
        name = y.pop(0)
        for idx, taxon in enumerate(y):
            y[idx] = taxon.split(':')[0]
        query_ham_dict[name] = y
    f.close()
    
    print ('{} seconds processing closest leaves'.format(time.perf_counter() - t0))

    #### 11.1.2023 - Chengze Shen
    # change q_dict to make sure its keys do not contain illegal characters
    # for directories
    pool = ProcessPoolExecutor(max_workers=num_threads, initializer=initializer,
            initargs=(tree, leaf_dict, ref_dict, query_ham_dict))

    qnum = 0
    for name in list(q_dict.keys()):
        seq = q_dict[name]
        q_dict.pop(name)
        q_dict[qnum] = (name, seq)
        qnum += 1

    # place each query sequence
    futures = []
    for qnum, val in q_dict.items():
        name, seq = val

        _args = (output, run, qnum, name, seq, n, info, subtree_flag, num_threads) 
        futures.append(pool.submit(runner, _args))

    for future in tqdm(concurrent.futures.as_completed(futures),
            total=len(futures), **tqdm_styles):
        tmp_output, time_extract_subtree, time_epa_ng = future.result()
        if tmp_output == None:
            continue

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
    shutil.rmtree("{}/tmp{}".format(args.outdir, run))
    
    
def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


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
                        help="Output file name", required=False, default="EPA-ng-SCAMPP")
    
    parser.add_argument("-m", "--model", type=str,
                        help="Model used for edge distances",
                        required=False, default="GTR")

    parser.add_argument("-b", "--subtreesize", type=int,
                        help="Integer size of the subtree",
                        required=False, default=2000)
    
    parser.add_argument("-s", "--subtreetype", type=str,
                        help="d (default) for edge weighted distances, n for node distances, h for hamming distances",
                        required=False, default="d")
    
    parser.add_argument("-n","--tmpfilenbr", type=int,
                        help="tmp file number",
                        required=False, default=0)
    
    parser.add_argument("-q", "--qalignment", type=str,
                        help="Path to query sequence alignment in fasta format (ref alignment separate)",
                        required=False, default="")
    
    parser.add_argument("-f", "--fragmentflag", type=str2bool,
                        help="boolean, True if queries contain fragments",
                        required=False, default=True)
    parser.add_argument('--threads', type=int,
                        help="Number of threads for EPA-ng",
                        required=False, default=1)

    parser.add_argument("-v", "--version", action="version", version="2.0.1", help="show the version number and exit")
                       
    return parser.parse_args()

if __name__ == "__main__":
    main(parseArgs())
