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

    # output path, ref, query, backbone tree, info

    aln_dict = utils.read_data(aln)
    ref_dict, q_dict = utils.seperate(aln_dict, query)
    tree = dendropy.Tree.get(path=tree_path, schema="newick", rooting='force-unrooted')
    
    tree.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    namespace = tree.taxon_namespace
    jplace = dict()
    placements = []

    counter = 0
    for e in tree.postorder_edge_iter():
        e.label = counter
        counter += 1
    jplace["tree"] = tree.as_string(schema="newick", suppress_edge_lengths=False, \
        edge_label_compose_fn=utils.edge_labeling, suppress_rooting=True)[:-1]

    files = []
    try: 
        os.mkdir("tmp{}".format(run))
    except:
    	print("tmp directory already exists")    
    try:
        os.mkdir(output)
    except:
    	print("path directory already exists")
        

    for name, seq in q_dict.items():
        tmp_tree = "../scratch/tmp{}/tree_".format(run)+name
        tmp_aln = "../scratch/tmp{}/aln".format(run)+name+".fa"
        tmp_output = output+"/"+name+".jplace"
        
        if not os.path.exists(tmp_output):
            y = utils.find_y(seq, ref_dict)
            nodes = utils.subtree_nodes(tree, y, n)
            subtree = tree.extract_tree_with_taxa(nodes)
        
            subtree.collapse_basal_bifurcation(set_as_unrooted_tree=True)
            subtree.write(path=tmp_tree, schema="newick", suppress_rooting=True)



        
      
            f = open(tmp_aln, "w")
            f.write(">"+name)
            f.write("\n")
            f.write(seq+"\n")
            for ns in nodes:
                f.write(">"+ns.label+"\n")
                f.write(ref_dict[ns.label])
                f.write("\n")

            f.close()

            os.system("./pplacer -m {} -s {} -t {} --keep-at-most 1 -o {} {}".format(model, info, tmp_tree, tmp_output, tmp_aln))


        place_file = open(tmp_output, 'r')
        place_json = json.load(place_file)

        if len(place_json["placements"]) > 0:

            added_tree = dendropy.Tree.get(data=place_json["tree"], 
            				    schema="newick", 
                                            is_parse_jplace_tokens = True,
        				    taxon_namespace=namespace)
        
            tmp_place = place_json["placements"][0]
            for i in range(len(tmp_place["p"])):
                edge_num = tmp_place["p"][i][1]
                edge_distal = tmp_place["p"][i][0]        
                

                place_edge = added_tree.edges(\
                        filter_fn=(lambda x: str(x.edge_number)==str(edge_num)))[0]
                right_n = place_edge.head_node
                left_n = right_n.parent_node
                
                #left and right path_l and path_r are in added_tree
                
                left, path_l = utils.find_closest(left_n, {left_n, right_n})
                
                
                right, path_r = utils.find_closest(right_n, {left_n, right_n})
                
                left = tree.find_node_with_taxon_label(left.taxon.label)
                right = tree.find_node_with_taxon_label(right.taxon.label)
                # now left right and path are in tree
                _, path = utils.find_closest(left, {left}, y=right)

                length = sum([x.length for x in path_l])+edge_distal# length through subtree
                target_edge = path[-1]
                
                for j in range(len(path)):
                    length -= path[j].length
                    if length < 0:
                        target_edge = path[j]
                        break
                                               
                tmp_place["p"][i][0] = 0

                tmp_place["p"][i][0] = target_edge.length+length
                tmp_place["p"][i][1] = target_edge.label
                                
            placements.append(tmp_place.copy())

        place_file.close()


    jplace["placements"] = placements
    jplace["metadata"] = {"invocation": " ".join(sys.argv)}
    jplace["version"] = 3
    jplace["fields"] = ["distal_length", "edge_num", "like_weight_ratio", \
            "likelihood", "pendant_length"]

    output = open(output+"/pplacer-XR.jplace", 'w')
    json.dump(jplace, output, sort_keys=True , indent=4)
    output.close()
    
    shutil.rmtree("tmp{}".format(run))
