__author__ = 'clipo'

import csv

import argparse

import logging as logger

import operator
import itertools

import os

import networkx as nx
import networkx.algorithms.isomorphism as iso
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from networkx.algorithms.isomorphism.isomorph import graph_could_be_isomorphic as isomorphic


class continuityAnalysis():

    def __init__(self):
        self.columnNames = []
        self.labels=[]
        self.taxa={}
        self.countOfTaxa=0
        self.graph=nx.Graph(name="MaximumGraph", GraphID=1, is_directed=False)
        self.minMaxGraph=nx.Graph(name="MaximumParsimony", GraphID=2, is_directed=False)
        self.outputDirectory =""
        self.FalseList=[None,0,False,"None","0","False"]

        self.dimensions={1:"Location of Maximum Blade Width",
            2: "Base Shape",
            3: "Basal-Indentation",
            4: "Constriction Ratio",
            5: "Outer Tang Angle",
            6: "Tang-Tip Shape",
            7: "Fluting",
            8: "Length/Width Ratio"}

        #print self.dimensions[2]

        self.classification={1:{1:"Proximal Quarter",2:"Secondmost Proximal Quarter",3:"Thirdmost Proximal Quarter",4:"Distal Quarter"},
		   2:{1:"Arc/Round",2:"Normal Curve",3:"Triangular", 4:"Folsomoid",5:"Flat",6:"Convex"},
		   3:{1:"No indentation",2:"Shallow",3:"Deep",4:"Very Deep"},
		   4:{1:"1.0",2:"0.90-0.99",3:"0.80-0.89",4:"0.70-0.79",5:"0.60-0.69",6:"0.50-0.59"},
		   5:{1:"93-115",2:"88-92",3:"81-87",4:"66-80",5:"51-65",6:"<=50"},
		   6:{1:"Pointed",2:"Round",3:"Blunt"},
		   7:{1:"Absent",2:"Present"},
		   8:{1:"1.00-1.99",2:"2.00-2.99",3:"3.00-3.99",4:"4.00-4.99",5:"5.00-5.99",6:">=6.00"}}

    def openFile(self, filename):
        try:
            logger.debug("trying to open: %s ", filename)
            file = open(filename, 'r')
        except csv.Error as e:
            logger.error("Cannot open %s. Error: %s", filename, e)
            sys.exit('file %s does not open: %s') % ( filename, e)

        reader = csv.reader(file, delimiter=' ', quotechar='|')
        for row in reader:
            row = map(str, row)
            label = row[0]
            self.labels.append(label)
            row.pop(0)
            self.taxa[label]=str(row[0])
            #print "characters: ", str(row[0])
            self.countOfTaxa += 1
        return True

    def saveGraph(self, graph, filename):
        nx.write_gml(graph, filename)

    def all_pairs(self, lst):
        return list((itertools.permutations(lst, 2)))

    def compareTaxa(self, taxa1, taxa2):
        #print self.taxa[taxa1], "-", self.taxa[taxa2]
        numberDiff = 0
        numberSame = 0
        count = 0
        dimensionsChanged =[]
        traitsChanged=[]
        traitsSame=[]
        dimensionsSame=[]
        for t in self.taxa[taxa1]:
            #print t, "-", self.taxa[taxa2][count]
            if t == self.taxa[taxa2][count]:
                numberSame += 1
                #print "dimensions: ",self.dimensions[count+1]
                dimensionsChanged.append(self.dimensions[count+1])
                change = self.dimensions[count+1]+":"+self.classification[count+1][int(self.taxa[taxa2][count])]+"->"+self.classification[count+1][int(t)]
                traitsChanged.append(change)
            else:
                #print "dimensions: ",self.dimensions[count+1]
                dimensionsChanged.append(self.dimensions[count+1])
                change = self.dimensions[count+1]+":"+self.classification[count+1][int(self.taxa[taxa2][count])]+"->"+self.classification[count+1][int(t)]
                traitsChanged.append(change)
            count += 1
        #print "1:  :", self.taxa[taxa1], "2: ", self.taxa[taxa2], "-->", number
        #print "dimensionsChanged: ", dimensionsChanged
        #print "traitsChanged: ", traitsChanged
        return  numberDiff, dimensionsChanged,traitsChanged,numberSame, dimensionsSame, traitsSame,

    def createGraph(self):
        allPairs = self.all_pairs(self.labels)
        for pairs in allPairs:
            num,dim,traits, numsame, dimsame, traitssame=self.compareTaxa(pairs[0],pairs[1])
            stuffChanged=str(dim)+"=>"+str(traits)
            #print "stuff:", stuffChanged
            if pairs[0] not in self.graph.nodes():
                self.graph.add_node(pairs[0], name=pairs[0], characterTraits=self.taxa[pairs[0]], connectedTo=pairs[1])
            if pairs[1] not in self.graph.nodes():
                self.graph.add_node(pairs[1], name=pairs[1], characterTraits=self.taxa[pairs[1]], connectedTo=pairs[0])

            self.graph.add_edge(pairs[0], pairs[1],
                                weight=self.compareTaxa(pairs[0],pairs[1]),
                                dims=dim,
                                traits=traits,
                                stuffChanged=stuffChanged)
            #print stuffChanged.strip('[]')

    def saveGraph(self, graph, filename):
        nx.write_gml(graph, filename[:-4]+".gml")

    ## from a "summed" graph, create a "min max" solution -- using Counts
    def createMinMaxGraph(self):
        ## first need to find the pairs with the maximum occurrence, then we work down from there until all of the
        ## nodes are included
        ## the weight

        maxWeight = 0
        pairsHash = {}
        traitList={}
        dimList={}
        stuffChanged={}
        for e in self.graph.edges_iter():
            d = self.graph.get_edge_data(*e)
            fromTaxa = e[0]
            toTaxa = e[1]
            #print d['weight']
            #print "weight: ", d['weight'][0]
            currentWeight = int(d['weight'][0])
            dimensions=d['weight'][1]
            traits=d['weight'][2]
            #stuff=d['weight'][3]
            pairsHash[fromTaxa + "*" + toTaxa] = currentWeight
            label = fromTaxa + "*" + toTaxa
            traitList[label]=traits
            dimList[label]=dimensions
            #stuffChanged[label]=stuff

        matchOnThisLevel = False
        currentValue = 0
        for key, value in sorted(pairsHash.iteritems(), key=operator.itemgetter(1), reverse=True):
            #print key, "->", value
            if value==0:
                value=.0000000000001
            if currentValue == 0:
                currentValue = value
            elif value < currentValue:
                matchOnThisLevel = False  ## we need to match all the connections with equivalent weights (otherwise we
                ## would stop after the nodes are included the first time which would be arbitrary)
                ## here we set the flag to false.
            taxa1, taxa2 = key.split("*")
            #print ass1, "-", ass2, "---",value
            if taxa1 not in self.minMaxGraph.nodes():
                self.minMaxGraph.add_node(taxa1, name=taxa1,characterTraits=self.taxa[taxa1])
            if taxa2 not in self.minMaxGraph.nodes():
                self.minMaxGraph.add_node(taxa2, name=taxa2, characterTraits=self.taxa[taxa2])
            if nx.has_path(self.minMaxGraph, taxa1, taxa2) == False or matchOnThisLevel == True:
                matchOnThisLevel = True   ## setting this true allows us to match the condition that at least one match was
                ## made at this level
                self.minMaxGraph.add_path([taxa1, taxa2], weight=value, dimensions=str(dimList[key]).strip('[]'),
                                          traits=str(traitList[key]).strip('[]'),
                                          #traitChanged=str(stuffChanged[key].strip('[]')),
                                          inverseweight=(1/value ))

        ## Output to file and to the screen
    def graphOutput(self):
        graph=self.minMaxGraph
        ## Now make the graphic for set of graphs
        plt.rcParams['text.usetex'] = False

        basefilename = os.path.basename(self.args['inputfile'])[:-4]


        newfilename = self.args['outputdirectory'] + "/" + basefilename + "-out.vna"
        gmlfilename = self.args['outputdirectory'] + "/" + basefilename + "-out.gml"
        self.saveGraph(graph, gmlfilename)
        plt.figure(newfilename, figsize=(8, 8))
        os.environ["PATH"] += ":/usr/local/bin:/usr/local/opt/graphViz/:"
        pos = nx.graphviz_layout(graph)
        edgewidth = []

        ### Note the weights here are biased where the *small* differences are the largest (since its max value - diff)
        weights = nx.get_edge_attributes(graph, 'weight')
        for w in weights:
            edgewidth.append(weights[w])
        maxValue = max(edgewidth)
        widths = []
        for w in edgewidth:
            widths.append(((maxValue - w) + 1) * 5)

        assemblageSizes = []
        sizes = nx.get_node_attributes(graph, 'size')
        #print sizes
        for s in sizes:
            #print sizes[s]
            assemblageSizes.append(sizes[s])
        nx.draw_networkx_edges(graph, pos, alpha=0.3, width=widths)
        sizes = nx.get_node_attributes(graph, 'size')
        nx.draw_networkx_nodes(graph, pos, node_color='w', alpha=0.4)
        nx.draw_networkx_edges(graph, pos, alpha=0.4, node_size=0, width=1, edge_color='k')
        nx.draw_networkx_labels(graph, pos, fontsize=10)
        font = {'fontname': 'Helvetica',
                'color': 'k',
                'fontweight': 'bold',
                'fontsize': 10}
        plt.axis('off')
        #plt.savefig(newfilename, dpi=150)
        plt.show()

    def checkMinimumRequirements(self):
        try:
            from networkx import graphviz_layout
        except ImportError:
            raise ImportError(
                "This function needs Graphviz and either PyGraphviz or Pydot. Please install GraphViz from http://www.graphviz.org/")
        if self.args['inputfile'] in self.FalseList:
            sys.exit("Inputfile is a required input value: --inputfile=../testdata/testdata.txt")

    def addOptions(self, oldargs):
        self.args = {'debug': None,  'inputfile': None, 'outputdirectory': None, 'pdf':None,'separator':None, 'missing':None, 'similarity':None,
                     'header':None}

        for a in oldargs:
            #print a
            self.args[a] = oldargs[a]

    def process(self,args):
        self.addOptions(args)
        self.checkMinimumRequirements()
        self.openFile(self.args['inputfile'])
        self.createGraph()
        self.createMinMaxGraph()
        self.graphOutput()
        #self.saveGraph(self.minMaxGraph,self.args['inputfile'])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Conduct a continuity analysis')
    #parser.add_argument('--debug', '-d', default=None, help='Sets the DEBUG flag for massive amounts of annotated output.')
    parser.add_argument('--inputfile','-f', required=True,
                        help="The file to be analyzed (.txt file) ")
    parser.add_argument('--outputdirectory', '-o', default=".", help="directory in which output files should be written.")
    #parser.add_argument('--separator','-s', default="tab",
                    #help="The type of separator between characters (space, tab, none) ")
    #parser.add_argument('--missing','-m',default=None, help='What to do with missing values (?) (e.g., estimate, none)')
    #parser.add_argument('--similarity','-si',default="similarity", help="Use similarity or dissimlarity")
    #parser.add_argument('--header','-hd', default=None, help='Whether or not there is a header (None, Yes)')

    args={}
    try:
        args = vars(parser.parse_args())
    except IOError, msg:
        parser.error(str(msg))
        sys.exit()

    ca = continuityAnalysis()

    results = ca.process(args)

''''
From the command line:

python ./continuityAnalysis.py --inputfile=../testdata/pfg.txt "


As a module:

from continuityAnalysis import continuityAnalysis

ca = continuityAnalysis()

args={}
args{'inputfile'}="../testdata/testdata-5.txt"
args{'screen'}=1
args{'debug'}=1
args('graphs'}=1

results = ca.process(args)

'''''