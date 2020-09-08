# Read Contents of File
with open('genome.txt', 'r') as file:
    genome = file.read()

# Make the string all uppercase for ease of use
genome = genome.upper()
cDNA = genome.upper()

# Replace the nucleotides for the complementary DNA strand
cDNA = cDNA.replace("A", "X")
cDNA = cDNA.replace("T", "A")
cDNA = cDNA.replace("X", "T")
cDNA = cDNA.replace("C", "X")
cDNA = cDNA.replace("G", "C")
cDNA = cDNA.replace("X", "G")

# Create DNA Object
DNA = (genome, cDNA)

# Blast primer #4
# ("Sequence, Starting Point, Ending point, GC Content")
fPrimer = ("GGTTTTGTCGTGCCTGGTTT", 298, 317, .5)
rPrimer = ("AGCAGCCAAAACACAAGCTG", 462, 443, .5) # Sequence is reversed 

# Print Sequence to replicate
print(DNA[0][fPrimer[1]:rPrimer[1]])
print(DNA[1][fPrimer[1]:rPrimer[1]])