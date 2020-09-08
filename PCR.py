with open('genome.txt', 'r') as file:
    genome = file.read()

cDNA = genome.upper()

cDNA = cDNA.replace("A", "X")
cDNA = cDNA.replace("T", "A")
cDNA = cDNA.replace("X", "T")
cDNA = cDNA.replace("C", "X")
cDNA = cDNA.replace("G", "C")
cDNA = cDNA.replace("X", "G")

#Blast primer #4
#("Sequence, Starting Point, Ending point, GC Content")
fPrimer = ("GGTTTTGTCGTGCCTGGTTT", 298, 317, .5)
rPrimer = ("AGCAGCCAAAACACAAGCTG", 462, 443, .5) #Sequence is reversed 

print(genome[fPrimer[1]:rPrimer[1]])
print(cDNA[fPrimer[1]:rPrimer[1]])