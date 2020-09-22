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
    This function will run a simulation of PCR using the submitted DNA segment and primers.  For each cycle, replicate
    each double stranded DNA tuple within the list.
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
    forward_compliment = find_compliment(forward_sequence)
    reverse_compliment = find_compliment(reverse_sequence)

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
                if reverse_compliment[1] in strand:
                    # start at the back -> front
                    # reverse strings to iterate through lists forwards for my mental health
                    strand_to_add = reverse_sequence[::-1]
                    reverse_strand = strand[::-1]
                    start_index = reverse_strand.find(reverse_compliment[1][::-1]) + len(reverse_sequence)
                    end_index = start_index + falloff_rate

                    # copy up to the length of the falloff_rate
                    for base in reverse_strand[start_index:end_index]:
                        strand_to_add = strand_to_add[:] + replaceDict[base]

                    # reverse string again for correct 5'-3' order
                    strand_to_add = strand_to_add[::-1]

                    # add to new strand to DNA pool
                    new_pair.append(strand_to_add)

                elif forward_compliment[1] in strand:
                    strand_to_add = forward_sequence[:]
                    start_index = strand.find(forward_compliment[1]) + len(forward_compliment[0])
                    end_index = start_index + falloff_rate

                    for base in strand[start_index:end_index]:
                        strand_to_add = strand_to_add + replaceDict[base]

                    # add to new strand to DNA pool
                    new_pair.append(strand_to_add)

                else:
                    new_pair.append(strand_to_add)

            new_pair = tuple(new_pair)
            dna_copied.append(new_pair)

        replicated_dna.extend(dna_copied)

    return replicated_dna

def find_statistics(replicated_dna):
    """
    Find statistics on replicated DNA. Finds strands, max strand length, min strand length, and average strand length.
    Finds average GC content of strands. Plots distributions of strands
    :param replicated_dna:
    :return:
    """

    segment_lengths = []
    gc_contents = []
    for pair in replicated_dna:
        for strand in pair:
            if strand != '':
                segment_lengths.append(len(strand))

                # Find GC contents of both strands
                num_of_c = strand.count('C')
                num_of_g = strand.count('G')
                gc_content = num_of_c + num_of_g
                gc_contents.append(gc_content)

    num_of_strands = len(segment_lengths)

    max_length = max(segment_lengths)
    min_length = min(segment_lengths)
    avg_length = sum(segment_lengths) / len(segment_lengths)
    avg_gc_content = sum(gc_contents) / len(gc_contents)

    hist = plt.hist(segment_lengths)
    plt.xlabel('Strand Lengths')
    plt.ylabel('Frequency')
    plt.title('Distribution of Strand Lengths')
    print('Total Strands found:', num_of_strands)
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

    replicated_DNA = run_PCR(DNA, fPrimer, rPrimer, cycles=20, fall_off_rate_base=180)
    find_statistics(replicated_DNA)
    print('PCR executed in: ', time.time() - start_time)
