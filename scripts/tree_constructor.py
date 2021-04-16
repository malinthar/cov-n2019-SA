from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio import Phylo
from Bio import Nexus
import matplotlib
import matplotlib.pyplot as plt
# import dendropy
# from dendropy.calculate import treecompare
from ete3 import Tree
# Before runing this class the sequence files(protein1.fasta, ..) should be aligned using clustalx.Use the output file location in drawTree method
    
def drawTrees():
    t1 = Phylo.read("./../data/alignment/ORF1ab_polyprotein/ORF1ab_polyprotein.ph", "newick")
    t2 = Phylo.read("./../data/alignment/surface_glycoprotein/surface_glycoprotein.ph", "newick")

    Phylo.draw(t1)
    Phylo.draw(t2)
    
def main():
    drawTrees()
 
if __name__ == "__main__":
    main()

