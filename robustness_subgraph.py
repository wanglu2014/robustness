#!/usr/bin/python

"""
This script performs robustness analysis on the given network, 
which involves removing nodes from the network at random, or in reverse 
order of centrality measures (degree, betweenness, closeness, and 
eigenvector), and comparing the size of the largest component in the 
network to the fraction of nodes removed.

Usage: python robustness.py <infile> <outfile> <recalculate>

where infile is the name of the network file in gml format, outfile is the 
name of the output (pdf) file in which the results of the analysis is 
saved, and recalculate (True of False) specifies if the targeted attack is 
simultaneous (False), or sequential (True).
"""

import igraph,networkx,numpy,operator,random,sys,pandas

def betweenness(infile, recalculate = False):
    """
    Performs robustness analysis based on betweenness centrality,  
    on the network specified by infile using sequential (recalculate = True) 
    or simultaneous (recalculate = False) approach. Returns a list 
    with fraction of nodes removed, a list with the corresponding sizes of 
    the largest component of the network, and the overall vulnerability 
    of the network.
    """

    #seperate network into sub
    g = networkx.read_gml(infile,label='id')
    fungigraph = g.subgraph([n for n,attrdict in g.node.items() if attrdict ['group']=='fungi'])
    bactigraph = g.subgraph([n for n,attrdict in g.node.items() if attrdict ['group'] == 'Bact' ] )
    lipidgraph = g.subgraph([n for n,attrdict in g.node.items() if attrdict ['group'] == 'lipid' ] )
    ##
    
    #remove node
    def rem_betw(subgraph,g):
        m = networkx.betweenness_centrality(subgraph)
        l = sorted(m.items(), key = operator.itemgetter(1), reverse = True)
        x = []
        y = []
        largest_component = max(networkx.connected_components(g), key = len)
        n = len(subgraph.nodes())
        x.append(0)
        y.append(len(largest_component) * 1. / n)
        R = 0.0
        for i in range(1, n):
            g.remove_node(l.pop(0)[0])
            if recalculate:#True, then restore all nodes. Therefore, False
                m = networkx.betweenness_centrality(g)
                l = sorted(m.items(), key = operator.itemgetter(1), 
                       reverse = True)
            largest_component = max(networkx.connected_components(g), key = len)
            #print(len(g.nodes()))
            x.append(i * 1. / n)
            R += len(largest_component) * 1. / n
            y.append(len(largest_component) * 1. / n)
        return x, y, 0.5 - R / n
    d = {"fungibet" : rem_betw(fungigraph,g),
    "bactibet": rem_betw(bactigraph,g),
    "lipidbet": rem_betw(lipidgraph,g)}
    return(pandas.DataFrame.from_dict(d,orient='index').transpose())

def degree(infile, recalculate = False):
    #seperate network into sub
    g = networkx.read_gml(infile,label='id')
    fungigraph = g.subgraph([n for n,attrdict in g.node.items() if attrdict ['group']=='fungi'])
    bactigraph = g.subgraph([n for n,attrdict in g.node.items() if attrdict ['group'] == 'Bact' ] )
    lipidgraph = g.subgraph([n for n,attrdict in g.node.items() if attrdict ['group'] == 'lipid' ] )
    ##
    #remove node
    def rem_degr(subgraph,g):
        m = networkx.degree_centrality(subgraph)
        l = sorted(m.items(), key = operator.itemgetter(1), reverse = True)
        x = []
        y = []
        largest_component = max(networkx.connected_components(g), key = len)
        n = len(subgraph.nodes())
        x.append(0)
        y.append(len(largest_component) * 1. / n)
        R = 0.0
        for i in range(1, n):
            g.remove_node(l.pop(0)[0])
            if recalculate:#True, then restore all nodes. Therefore, False
                m = networkx.betweenness_centrality(g)
                l = sorted(m.items(), key = operator.itemgetter(1), 
                       reverse = True)
            largest_component = max(networkx.connected_components(g), key = len)
            #print(len(g.nodes()))
            x.append(i * 1. / n)
            R += len(largest_component) * 1. / n
            y.append(len(largest_component) * 1. / n)
        return x, y, 0.5 - R / n
    d = {"fungideg" : rem_degr(fungigraph,g),
    "bactideg": rem_degr(bactigraph,g),
    "lipiddeg": rem_degr(lipidgraph,g)}
    from IPython.core.debugger import Tracer; Tracer()();set_trace()
    return(pandas.DataFrame.from_dict(d,orient='index').transpose())

def rand(infile):
    #seperate network into sub
    g = networkx.read_gml(infile,label='id')
    fungigraph = g.subgraph([n for n,attrdict in g.node.items() if attrdict ['group']=='fungi'])
    bactigraph = g.subgraph([n for n,attrdict in g.node.items() if attrdict ['group'] == 'Bact' ] )
    lipidgraph = g.subgraph([n for n,attrdict in g.node.items() if attrdict ['group'] == 'lipid' ] )
    ##
        #remove node
    def rem_rand(subgraph,g):
        l = [(node, 0) for node in subgraph.nodes()]
        random.shuffle(l)
        x = []
        y = []
        largest_component = max(networkx.connected_components(g), key = len)
        n = len(subgraph.nodes())
        x.append(0)
        y.append(len(largest_component) * 1. / n)
        R = 0.0
        for i in range(1, n):
            g.remove_node(l.pop(0)[0])
            largest_component = max(networkx.connected_components(g), key = len)
            x.append(i * 1. / n)
            R += len(largest_component) * 1. / n
            y.append(len(largest_component) * 1. / n)
        return x, y, 0.5 - R / n
    d = {"fungig" : rem_rand(fungigraph,g),
    "bactig": rem_rand(bactigraph,g),
    "lipidg": rem_rand(lipidgraph,g)}
    return(pd.DataFrame.from_dict(d,orient='index').transpose())



def main(argv):
    """
    Entry point.
    """

    if len(argv) != 3:
        print "python robustness.py <infile> <outfile> <recalculate>"
        sys.exit(0)

    infile = argv[0]
    outfile = argv[1]
    if argv[2] == "True":
        recalculate = True
    else:
        recalculate = False
    VD = degree(infile, recalculate)
    VB = betweenness(infile, recalculate)
    #from IPython.core.debugger import Tracer; Tracer()();set_trace()
    #VR = rand(infile)
    

    df_c = pandas.concat([VD, VB], axis=1)

    df_c.to_csv(outfile)
if __name__ == "__main__":
    main(sys.argv[1:])