def find_compliment(genome):
    """
    Find the compliment to submitted RNA genome
    :param genome: String of RNA
    :return: Tuple of original string and complement
    """
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
    return DNA


def run_PCR(dna, forward_primer, reverse_primer, cycles=10):
    """
    This function will run a simulation of PCR using the submitted DNA segment and primers. First, input DNA will be
    will be shortened to only the section that will be replicated and saved to a tuple that is appended to an empty list.
    Loop through the number of cycles. For each cycle, replicate each double stranded DNA tuple within the list.
    Within each double stranded tuple, replicate each strand. Strand replication is done by finding which primer is the
    compliment to the input strand. The rest of the strand is then looped through, finding the rest of the base pairs.
    One large list of DNA tuples is then exported as a result.
    :param dna: Tuple with original RNA strand and cDNA strand
    :param forward_primer: Tuple for forward primer. Contains: sequence, starting point, ending point, GC content
    :param reverse_primer: Tuple for reverse primer. Contains: sequence, starting point, ending point, GC content
    :param cycles: Number of cycles to simulate. Defaults to 10.
    :return: list with each entry being half of a replicated DNA string. Next entry is the other half of the strand
    """

    forward_sequence = forward_primer[0]
    forward_start = forward_primer[1]
    reverse_sequence = reverse_primer[0][::-1]
    reverse_start = reverse_primer[1]

    # Shorten input DNA to focus only on the section we care about
    shortened_dna = list()
    shortened_dna.append(dna[0][forward_start:reverse_start])
    shortened_dna.append(dna[1][forward_start:reverse_start])

    # Make copy of DNA to store the results of the PCR
    replicated_dna = list()
    replicated_dna.append(tuple(shortened_dna))

    # Replicate the section
    for cycle in range(1, (cycles + 1)):
        dna_copied = []  # All of the new DNA pairs that will be found in the cycle
        for strands in replicated_dna:
            new_pair = list()  # New double strand of DNA (tuple)
            for strand in strands:  # For loop definition does the equivalent of denaturation
                strand_to_add = ''
                if forward_sequence[1:] in strand:
                    # start at the back -> front
                    # reverse strings to iterate through lists forwards for my mental health
                    strand_to_add = reverse_sequence[::-1]
                    reverse_strand = strand[::-1]

                    # for base in replicated_dna[0][-len(reverse_sequence) - 1::-1]: (this loops through it backwards)
                    for base in reverse_strand[len(reverse_sequence):]:
                        if base == "A":
                            strand_to_add = strand_to_add[:] + "T"
                        if base == "T":
                            strand_to_add = strand_to_add[:] + "A"
                        if base == "G":
                            strand_to_add = strand_to_add[:] + "C"
                        if base == "C":
                            strand_to_add = strand_to_add[:] + "G"

                    # reverse string again for correct 5'-3' order
                    strand_to_add = strand_to_add[::-1]

                    # add to new strand to DNA pool
                    new_pair.append(strand_to_add)

                elif reverse_sequence in strand:

                    strand_to_add = forward_sequence[1:]

                    # for base in replicated_dna[0][-len(reverse_sequence) - 1::-1]: (this loops through it backwards)
                    for base in strand[len(forward_sequence[1:]):]:
                        if base == "A":
                            strand_to_add = strand_to_add + "T"
                        if base == "T":
                            strand_to_add = strand_to_add + "A"
                        if base == "G":
                            strand_to_add = strand_to_add + "C"
                        if base == "C":
                            strand_to_add = strand_to_add + "G"

                    # add to new strand to DNA pool
                    new_pair.append(strand_to_add)

                else:
                    print('Neither primer was found in the strand. Something went wrong')

            new_pair = tuple(new_pair)
            dna_copied.append(new_pair)

        replicated_dna.extend(dna_copied)

        '''
        print("DNA after PCR")
        for rna in replicated_dna:
            print(rna)
        '''

    return replicated_dna


if __name__ == '__main__':
    # Read Contents of File
    with open('genome.txt', 'r') as file:
        genome = file.read()

    # Make the string all uppercase for ease of use
    genome = genome.upper()

    DNA = find_compliment(genome)

    # Blast primer #4
    # ("Sequence, Starting Point, Ending point, GC Content")
    fPrimer = ("GGTTTTGTCGTGCCTGGTTT", 298, 317, .5)
    rPrimer = ("AGCAGCCAAAACACAAGCTG", 462, 443, .5)  # Sequence is reversed

    # Print Sequence to replicate
    print(DNA[0][fPrimer[1]:rPrimer[1]])
    print(DNA[1][fPrimer[1]:rPrimer[1]])

    replicated_DNA = run_PCR(DNA, fPrimer, rPrimer, cycles=10)
