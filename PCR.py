import time
import random
import matplotlib.pyplot as plt
import matplotlib
import tkinter
matplotlib.use('tkagg')

# Find the Compliment
def find_compliment(RNA, flag):
    """
    Find the compliment to submitted RNA
    :param RNA: String of RNA
    :param flag: determine how to return, As a string, or a Tuple
    :return: Tuple of original string and complement, or a 
             single compliment string
    """
    cDNA = RNA.upper()
    # Replace the nucleotides for the complementary DNA strand
    cDNA = cDNA.replace("A", "X")
    cDNA = cDNA.replace("T", "A")
    cDNA = cDNA.replace("X", "T")
    cDNA = cDNA.replace("C", "X")
    cDNA = cDNA.replace("G", "C")
    cDNA = cDNA.replace("X", "G")

    if flag == 1:
        # Return Single String
        return cDNA
    else:
        # Create DNA Object
        DNA = (RNA, cDNA)
        return DNA


def run_PCR(dna, forward_primer, reverse_primer, cycles=10, falloff_base=180):
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
    :param falloff_base: The base falloff rate to be incremented using a random int between -50 and 50.
    :return: list with each entry being half of a replicated DNA string. Next entry is the other half of the strand
"""

    replaceDict = {"A": "T", "T": "A", "G": "C", "C": "G"}
    forward_sequence = forward_primer[0]
    reverse_sequence = reverse_primer[0][::-1]
    forward_compliment = find_compliment(forward_sequence, 0)
    reverse_compliment = find_compliment(reverse_sequence, 0)
    # Make copy of DNA to store the results of the PCR
    replicated_dna = [dna]
    # Replicate the section
    for cycle in range(1, (cycles + 1)):
        start_time = time.time()
        print(f"Cycle: {cycle}")
        # All of the new DNA pairs that will be found in the cycle
        dna_copied = list() 
        # Start Of Denaturation Step: Separate DNA into 2 Strands
        for dna in replicated_dna:
            # New double strand of DNA (convert to tuple at end of iteration)
            new_pair = list()  
            # For loop definition does the equivalent of denaturation
            for strand in dna:  
                # Calculate fall off rate
                falloff_rate = falloff_base + random.randint(-50,50)
                # Start of Annealing Step: Add the correct primer to the strand
                # Search for the reverse primer compliment in the strand
                if reverse_compliment[1] in strand:
                    # Reverse strings to iterate through lists forward -> back
                    # Add the reverse primer to a str variable 
                    strand_to_add = reverse_sequence[::-1]
                    
                    reverse_strand = strand[::-1]
                    start_index = reverse_strand.find(
                        reverse_compliment[1][::-1]) + len(reverse_sequence)
                    end_index = start_index + falloff_rate

                    # Start of Elongation Step: Attach
                    # copy up to the length of the falloff_rate
                    taq_polymerase = find_compliment(reverse_strand[start_index:end_index], 1)
                    strand_to_add = strand_to_add + taq_polymerase
                    # reverse string again for correct 5'-3' order
                    strand_to_add = strand_to_add[::-1]
                    # add to new strand to DNA pool
                    new_pair.append(strand_to_add)

                # Otherwise search for the forward primer compliment in the strand
                elif forward_compliment[1] in strand:
                    # Add the forward primer to a str variable
                    strand_to_add = forward_sequence[:]

                    start_index = strand.find(forward_compliment[1]) + len(forward_compliment[0])
                    end_index = start_index + falloff_rate

                    # Start of Elongation Step: Attach
                    # copy up to the length of the falloff_rate
                    taq_polymerase = find_compliment(strand[start_index:end_index], 1)
                    strand_to_add = strand_to_add + taq_polymerase
                    # add to new strand to DNA pool
                    new_pair.append(strand_to_add)

                # If not found then the primers have fallen off so dont replicate this strand
                else:
                    new_pair.append('')

            new_pair = tuple(new_pair)
            dna_copied.append(new_pair)

        replicated_dna.extend(dna_copied)
        print(f"Time To Complete: {time.time() - start_time}")
    return replicated_dna


def find_statistics(replicated_dna):
    """
    Find statistics on replicated DNA. Finds strands, 
    max strand length, min strand length, and average strand length.
    Finds average GC content of strands. Plots distributions of strands
    :param replicated_dna:
    :return:
    """
    segment_lengths = []
    gc_contents = []
    for pair in replicated_dna[1:]:
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
    avg_gc_content = (sum(gc_contents) / len(gc_contents)) / avg_length
    
    hist = plt.hist(segment_lengths)
    plt.xlabel('Strand Lengths')
    plt.ylabel('Frequency')
    plt.title('Distribution of Strand Lengths')
    print(f'Total Strands found: {num_of_strands}')
    print(f'Average GC Content: {"%0.2f" % (avg_gc_content*100)}%', )
    print(f'Max Length: {max_length}')
    print(f'Min Length: {min_length}')
    print(f'Average Length: {avg_length}')
    plt.show()
    return


if __name__ == '__main__':
    random.seed(99)
    start_time = time.time()

    # Setup Step 1: Read Contents of File
    with open('genome.txt', 'r') as file:
        RNA = file.read()

    # Make the string all uppercase for ease of use
    RNA = RNA.upper()
    # Setup Step 2: Find Compliments
    DNA = find_compliment(RNA, 0)

    # Setup Step 3: Define Primers
    # ("Sequence, Starting Point, Ending point, GC Content")

    fPrimer = ("GGTTTTGTCGTGCCTGGTTT", 297, 317, .5)
    rPrimer = ("AGCAGCCAAAACACAAGCTG", 464, 443, .5)  # Sequence is reversed

    replicated_DNA = run_PCR(DNA, fPrimer, rPrimer,
                             cycles=20, falloff_base=180)
    print('PCR executed in: ', time.time() - start_time)
    find_statistics(replicated_DNA)
