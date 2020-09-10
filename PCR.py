def run_PCR(dna, forward_primer, reverse_primer, cycles=10):
    """
    This function will run a simulation of PCR using the submitted DNA segment and primers.
    :param dna: Tuple with original RNA strand and cDNA strand
    :param forward_primer: Tuple for forward primer. Contains: sequence, starting point, ending point, GC content
    :param reverse_primer: Tuple for reverse primer. Contains: sequence, starting point, ending point, GC content
    :param cycles: Number of cycles to simulate. Defaults to 10.
    :return: list with each entry being half of a replicated DNA string. Next entry is the other half of the strand
    """

    forward_sequence = forward_primer[0]
    forward_start = forward_primer[1]
    forward_end = forward_primer[2]
    reverse_sequence =  reverse_primer[0][::-1]
    reverse_start = reverse_primer[1]
    reverse_end = reverse_primer[2]

    # Make copy of DNA to store the results of the PCR
    replicated_dna = list(dna)

    # Find only section that will be replicated
    replicated_dna[0] = replicated_dna[0][forward_start:reverse_start]
    replicated_dna[1] = replicated_dna[1][forward_start:reverse_start]

    # Replicate the section
    for cycle in range(1, (cycles + 1)):
        dna_copied = []
        for i, strand in zip(range(0, (len(replicated_dna) + 1)), replicated_dna):
            strand_to_add = ''
            if forward_sequence in strand:
                # TODO: Loop through strand backwards, checking and adding opposite bases
                strand_to_add = reverse_sequence

            elif reverse_sequence in strand:
                # TODO: Loop through strand forwards, checking and adding opposite bases
            else:
                print('Neither primer was found in the strand. Something went wrong')

    return replicated_dna

if __name__ == '__main__':
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

    '''
    #replicated_DNA = run_PCR(DNA, fPrimer, rPrimer, cycles=1)

    for entry in replicated_DNA:
        print(entry)
    '''