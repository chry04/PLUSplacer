from dendropy import *
import numpy as np
import heapq



def read_data(aln):
    """ Load the query and reference sequence from the alignment file
    
    Parameters
    ----------
    aln : multiple sequence alignment containing reference taxa and query sequence
    
    Returns
    -------
    dictionary containing sequences with taxon label keys
    
    """
    
    f = open(aln)
    result = dict()

    taxa = ""
    seq = ""
    for line in f:
        if line[0] == '>':
            if taxa != "":
                result[taxa] = seq
            taxa = line[1:-1]
            seq = ""

        elif line == "/n":
            continue
        else:
            seq += line[:-1]
            
    if taxa != "":
        result[taxa] = seq


    return result

def seperate(aln_dict, query):
    """ Separate the query sequences from the reference sequences
    
    Parameters
    ----------
    aln_dict : Sequence dictionary with taxon label keys
    
    Returns
    -------
    separate dictionaries containing query sequences and referece sequences with taxon label keys
    
    """
    
    f = open(query)
    q_name = set()
    for line in f:
        q_name.add(line[:-1])

    ref = dict()
    query = dict()

    for key, value in aln_dict.items():
        if key in q_name:
            query[key] = value
        else:
            ref[key] = value
    
    return ref, query

def hamming(seq1, seq2):
    """ Returns hamming distance between two sequences
    
    Parameters
    ----------
    seq1 : query sequence
    seq2 : reference sequence
    
    Returns
    -------
    integer hamming distance between query sequence and reference sequence
    
    """
    return len([1 for i in range(len(seq1)) if seq1[i] != seq2[i]])

def process_backbone_tree(aln, output_folder, n_ref, n_query, n_set):
    """ Not used ~ Samples alignment for dataset creation
    and creates fasta file where sites containing 95% gaps are removed.
    
    output n dataset with
    each dataset have a aln(alignment for both query and reference),
    query.txt(taxon name for query)
    """
    aln_dict = read_data(aln)
    n_row = len(aln_dict)
    n_col = len(aln_dict.key[0])
    array = np.zeros((n_row, n_col))
    taxons = []
    counter = 0

    for key, value in aln_dict.items():
        taxons.append(key)
        array[counter,:] = value

    for i in range(n_set):
        choices = np.random.choice(n_row, size=(n_ref+n_query), replace=False)
        query = choices[:n_query]
        aln_matrix = array[choices]

        kept = []
        for i in range(n_col):
            dashes = float((np.char.count(aln_matrix[:, i], '-')[0])/n_row)

            if dashes < 0.95:
                kept.append(i)

        aln = aln_martix[:,kept]

        aln_path = output+"/aln"+i+".fa"
        q_path = output+"/query"+i

        f = open(aln_path, "w")
        for j in range(len(choices)):
            f.write(">"+taxons[counter[i]]+"\n")
            f.write(aln_martix[i])
            f.write("\n")
        f.close()

        f = open(q_path, "w")
        for j in range(len(query)):
            f.write(taxons[query[i]])
            f.write("\n")
        f.close()
        

def find_y(x, ref):
     """ Returns leaf label for closest sister taxon l
    
    Parameters
    ----------
    x : aligned query sequence
    ref : reference multiple sequence alignment dictionary 
    
    Returns
    -------
    leaf label for taxon with smallest hamming distacne to query sequence
    
    """
    low = len(x)
    y = ""
    for name, seq in ref.items():
        h_dist = hamming(x, seq)
        if h_dist < low:
            low = h_dist
            y = name
    return y


def find_closest(x, visited, y=None):
    """ Returns leaf label for closest leaf to the node x through path not travelling through visited.
    If y is populated returns path from x to y not travelling through nodes in visited.
    
    Parameters
    ----------
    x : dendropy node object
    visited : list containing dendropy node objects 
    y : dendropy node object
    
    Returns
    -------
    If y == None : dendropy node object of closest leaf y to the node x through path not travelling through nodes in visited, 
                   list containing dendropy node objects on path to that leaf y from node x
    If y != None : dendropy node object y, 
                   list containing dendropy node objects on path from node x to leaf y not travelling through nodes in visited
    
    """
    queue = []
    counter = 1
    visited.add(x)

    if x.parent_node and x.parent_node not in visited:
        tmp = []
        tmp.append(x.edge)
        heapq.heappush(queue, [x.edge_length, counter,\
            tmp, x.parent_node])
        counter += 1


    for child in x.child_nodes():
        if child and child not in visited:
            tmp = []
            tmp.append(child.edge)
            heapq.heappush(queue, [child.edge_length, counter,\
                    tmp, child])
            counter += 1

    while len(queue) > 0:
        try:
            [length, _, path, node] = heapq.heappop(queue)
        except IndexError:
            break

        visited.add(node)
        if node.is_leaf():
            if (not y) or node.taxon.label==y.taxon.label:
                return node, path
            else: 
                continue

                    
        if node.parent_node and node.parent_node not in visited:
            tmp = path.copy()
            tmp.append(node.edge)
            heapq.heappush(queue, [length+node.edge_length, counter,\
                tmp, node.parent_node])
            counter += 1


        for child in node.child_nodes():
            if child and child not in visited:
                tmp = path.copy()
                tmp.append(child.edge)
                heapq.heappush(queue, [length+child.edge_length, counter,\
                        tmp, child])
                counter += 1

    return x, [x.edge]


def subtree_nodes(tree, y, n):
    """ Returns list of length n of taxon dendropy node objects (leaves) closest to sister taxon
    
    Parameters
    ----------
    tree : dendropy tree object
    y : taxon label for closest sister taxon
    n = number of taxa contained in subtree
    
    Returns
    -------
    list of dendropy node objects corresponding to leaves in the subtree
    
    """
    leaf_y = tree.find_node_with_taxon_label(y)
    queue = [(0, 0, leaf_y.parent_node)]
    
    leaves = [leaf_y]
    visited = {leaf_y}

    counter = 1

    while len(leaves) < n:
        try:
            (length, _, node) = heapq.heappop(queue)
        except IndexError:
            break

        visited.add(node)
        if node.is_leaf():
            leaves.append(node)

        for child in node.adjacent_nodes():
            if child not in visited:
                heapq.heappush(queue, (length+1, counter, child))
                counter += 1

    
    result = []
    for item in leaves:
        result.append(item.taxon)

    return result

def subtree_nodes_with_edge_length(tree, y, n):
    """ Not used ~ Returns list of length n of taxon dendropy node objects (leaves) closest to sister taxon
    
    Parameters
    ----------
    tree : dendropy tree object
    y : taxon label for closest sister taxon
    n = number of taxa contained in subtree
    
    Returns
    -------
    list of dendropy node objects corresponding to leaves in the subtree
    
    """

    leaf_y = tree.find_node_with_taxon_label(y)
    queue = [(leaf_y.edge_length, 0, leaf_y.parent_node)]
    
    leaves = [leaf_y]
    visited = {leaf_y}

    counter = 1

    while len(leaves) < n:
        try:
            (length, _, node) = heapq.heappop(queue)
        except IndexError:
            break

        visited.add(node)
        if node.is_leaf():
            leaves.append(node)

        for child in node.adjacent_nodes():
            if child not in visited:
                heapq.heappush(queue, (length+child.edge_length, counter, child))
                counter += 1

    
    result = []
    for item in leaves:
        result.append(item.taxon)

    return result

def edge_labeling(l):
    """ Returns jplace tokens attached to edge lengths
      
    Parameters
    ----------
    l : dendropy edge object
    
    Returns
    -------
    string containing jplace tokens attached to edge lengths
    """
    return str(l.length)+"{"+str(l.label)+"}"



def compareTreesFromPath(treePath1, treePath2):
    """
    This code was written by Erin Molloy and was borrowed from :
    https://github.com/ekmolloy/njmerge/blob/master/python/compare_trees.py
    """
    print("Comparing {} with {}".format(treePath1, treePath2))

    tax = TaxonNamespace()
    tr1 = Tree.get(path=treePath1,
                            schema='newick',
                            rooting='force-unrooted',
                            taxon_namespace=tax,
                            preserve_underscores=True)
    tr2 = Tree.get(path=treePath2,
                            schema='newick',
                            rooting='force-unrooted',
                            taxon_namespace=tax,
                            preserve_underscores=True)

    tr1.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    tr2.collapse_basal_bifurcation(set_as_unrooted_tree=True)

    return compareDendropyTrees(tr1, tr2)
    # print("RF distance on %d shared leaves: %d" % (nl, fp + fn))


def compareDendropyTrees(tr1, tr2):
    """
    This code was written by Erin Molloy and was borrowed from :
    https://github.com/ekmolloy/njmerge/blob/master/python/compare_trees.py
    """
    from dendropy.calculate.treecompare \
        import false_positives_and_negatives

    lb1 = set([l.taxon.label for l in tr1.leaf_nodes()])
    lb2 = set([l.taxon.label for l in tr2.leaf_nodes()])

    com = lb1.intersection(lb2)
    if com != lb1 or com != lb2:
        com = list(com)
        tns = TaxonNamespace(com)

        tr1.retain_taxa_with_labels(com)
        tr1.migrate_taxon_namespace(tns)

        tr2.retain_taxa_with_labels(com)
        tr2.migrate_taxon_namespace(tns)
    com = list(com)

    tr1.update_bipartitions()
    tr2.update_bipartitions()

    nl = len(com)
    ei1 = len(tr1.internal_edges(exclude_seed_edge=True))
    ei2 = len(tr2.internal_edges(exclude_seed_edge=True))

    [fp, fn] = false_positives_and_negatives(tr1, tr2)
    rf = float(fp + fn) / (ei1 + ei2)

    return (nl, ei1, ei2, fp, fn, rf)

def writePhylip(alignment, filePath, taxa = None):
    """  writes file in Phylip format
    This code was provided by Vladimir Smirnov in CS 581 at UIUC 
      
    Parameters
    ----------
    alignment : alignment dictionary where key is taxon name and value is object with "seq" attribute
    filePath : path for output Phylip file 
    taxa : taxa to be included
    
    """
    
    maxChars = 0
    lines = []
    for tag in alignment:
        if taxa is None or tag in taxa:
            lines.append("{} {}\n".format(tag, alignment[tag]))
            maxChars = max(maxChars, len(alignment[tag]))
    with open(filePath, 'w') as textFile:
        textFile.write("{} {}\n".format(len(lines), maxChars))
        for line in lines:
            textFile.write(line)
