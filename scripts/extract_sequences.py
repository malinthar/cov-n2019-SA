import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# read proteins table excel file. Provide he correct location
xls = pd.read_excel('./../data/protein_table.xlsx')

# The list of proteins that we consider
protein_set = ["ORF1ab_polyprotein", "surface_glycoprotein"]

protein1 = []
protein2 = []

# accession numbers of each common species.
acc_numbers = ["MT371047.1","MT371048.1","MT371049.1","MT371050.1","MT396246.2","MW173254.1","MW242999.1",
"MW645475.1","MT262993.1","MW447627.1","MW411961.1","MW242667.1","MT502774.1","MT913012.1",
"MW532098.1","MW532101.1"]

positions_ORF1ab_polyprotein = {}
positions_surface_glycoprotein = {}

def findPositions():
    print(xls)
    for i in range (0,2):
        for j in range(0,16):
            if(i==0):
                pos=[xls["ORF1ab_polyprotein_start"][j],xls["ORF1ab_polyprotein_end"][j]]
                positions_ORF1ab_polyprotein[acc_numbers[j]] = pos
            else:
                pos=[xls["surface_glycoprotein_start"][j],xls["surface_glycoprotein_end"][j]]
                positions_surface_glycoprotein[acc_numbers[j]] = pos


def extractGeneSequences():
    for acc in acc_numbers:
        fileName = "./../data/genomes/all/"+acc+".fasta"
        with open(fileName) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                for i in range (0,2):
                    if(i == 0):
                        start = positions_ORF1ab_polyprotein[acc][0]
                        stop = positions_ORF1ab_polyprotein[acc][1]
                        temp = [record.id, record.seq[start:stop]]
                        protein1.append(temp)
                    else:
                        start = positions_surface_glycoprotein[acc][0]
                        stop = positions_surface_glycoprotein[acc][1]
                        temp = [record.id, record.seq[start:stop]]
                        protein2.append(temp)


def writeSequences():
    for i in range(0, 2):
        fileName = "./../data/output/"+ protein_set[i]+".fasta"
        if(i == 0):
            seqlist = protein1
            
        else:
            seqlist = protein2

        sequences = []
        for seq in seqlist:
            record = SeqRecord(seq[1], id=seq[0], name=seq[0],
                               description=protein_set[i]+" protein of "+seq[0])
            sequences.append(record)
        SeqIO.write(sequences, fileName, "fasta")


def main():
    findPositions()
    print(positions_ORF1ab_polyprotein)
    print(positions_surface_glycoprotein)
    extractGeneSequences()
    writeSequences()
    print("Successfully extracted sequeces and saved in ./data/output directory")

if __name__ == "__main__":
    main()
