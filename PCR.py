import time
import random
import matplotlib.pyplot as plt

def find_compliment(RNA):
    """
    Find the compliment to submitted RNA RNA
    :param RNA: String of RNA
    :return: Tuple of original string and complement
    """
    cDNA = RNA.upper()

    # Replace the nucleotides for the complementary DNA strand
    cDNA = cDNA.replace("A", "X")
    cDNA = cDNA.replace("T", "A")
    cDNA = cDNA.replace("X", "T")
    cDNA = cDNA.replace("C", "X")
    cDNA = cDNA.replace("G", "C")
    cDNA = cDNA.replace("X", "G")

    # Create DNA Object
    DNA = (RNA, cDNA)
    return DNA


def run_PCR(dna, forward_primer, reverse_primer, cycles=10, fall_off_rate_base=180):
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
    :param fall_off_rate_base: The base falloff rate to be incremented using a random int between -50 and 50.
    :return: list with each entry being half of a replicated DNA string. Next entry is the other half of the strand
    """

    replaceDict = {"A": "T", "T": "A", "G": "C", "C": "G"}
    forward_sequence = forward_primer[0]
    reverse_sequence = reverse_primer[0][::-1]

    # Make copy of DNA to store the results of the PCR
    replicated_dna = list()
    replicated_dna.append(tuple(dna))

    # Replicate the section
    for cycle in range(1, (cycles + 1)):
        dna_copied = []  # All of the new DNA pairs that will be found in the cycle
        for strands in replicated_dna:
            new_pair = list()  # New double strand of DNA (convert to tuple at end of iteration)
            for strand in strands:  # For loop definition does the equivalent of denaturation
                strand_to_add = ''
                falloff_rate = fall_off_rate_base + random.randint(-50, 50)
                print(falloff_rate)
                #if forward_sequence[1:] in strand:
                if forward_sequence in strand:
                    # start at the back -> front
                    # reverse strings to iterate through lists forwards for my mental health
                    strand_to_add = reverse_sequence[::-1]
                    reverse_strand = strand[::-1]

                    # copy up to the length of the falloff_rate
                    for base in reverse_strand[len(reverse_sequence):]:
                        if ((len(strand_to_add) - 20) <= falloff_rate):
                            strand_to_add = strand_to_add[:] + replaceDict[base]

                    # reverse string again for correct 5'-3' order
                    strand_to_add = strand_to_add[::-1]

                    # add to new strand to DNA pool
                    new_pair.append(strand_to_add)

                elif reverse_sequence in strand:
                    strand_to_add = forward_sequence[:]

                    for base in strand[len(forward_sequence[:]):]:
                        if ((len(strand_to_add) - 20) <= falloff_rate):
                            strand_to_add = strand_to_add + replaceDict[base]

                    # add to new strand to DNA pool
                    new_pair.append(strand_to_add)

                else:
                    # TODO: Put in flagging mechanism
                    print('Neither primer was found in the strand. Something went wrong')

            new_pair = tuple(new_pair)
            dna_copied.append(new_pair)

        replicated_dna.extend(dna_copied)
        print(f"Cycle: {cycle}")

        '''
        print("DNA after PCR")
        for rna in replicated_dna[1:]:
            print(rna)
            '''


    return replicated_dna

def find_statistics(replicated_dna):

    segment_lengths = []
    gc_contents = []
    for pair in replicated_dna:
        segment_lengths.append(len(pair[0]))
        segment_lengths.append(len(pair[1]))

        # Find GC contents of both strands
        num_of_c = pair[0].count('C')
        num_of_g = pair[0].count('G')
        gc_content = num_of_c + num_of_g
        gc_contents.append(gc_content)
        num_of_c = pair[1].count('C')
        num_of_g = pair[1].count('G')
        gc_content = num_of_c + num_of_g
        gc_contents.append(gc_content)

    num_of_strands = len(replicated_dna) * 2

    max_length = max(segment_lengths)
    min_length = min(segment_lengths)
    avg_length = sum(segment_lengths) / len(segment_lengths)
    avg_gc_content = sum(gc_contents) / len(gc_contents)

    hist = plt.hist(segment_lengths)
    plt.xlabel('Strand Lengths')
    plt.ylabel('Frequency')
    plt.title('Distribution of Strand Lengths')
    print('Average GC Content:', avg_gc_content)
    print('Max Length:', max_length)
    print('Min Length:', min_length)
    print('Average Length:', avg_length)
    plt.show()

    return

if __name__ == '__main__':
    #random.seed(99)
    start_time = time.time()

    # Read Contents of File
    with open('genome.txt', 'r') as file:
        RNA = file.read()

    # Make the string all uppercase for ease of use
    RNA = RNA.upper()

    DNA = find_compliment(RNA)

    # Blast primer #4
    # ("Sequence, Starting Point, Ending point, GC Content")
    fPrimer = ("GGTTTTGTCGTGCCTGGTTT", 297, 317, .5)
    rPrimer = ("AGCAGCCAAAACACAAGCTG", 464, 443, .5)  # Sequence is reversed

    # Print Sequence to replicate
    #print(DNA[0][fPrimer[1]:rPrimer[1]])
    #print(DNA[1][fPrimer[1]:rPrimer[1]])

    replicated_DNA = run_PCR(DNA, fPrimer, rPrimer, cycles=15, fall_off_rate_base=180)
    find_statistics(replicated_DNA)
    print('PCR executed in: ', time.time() - start_time)

