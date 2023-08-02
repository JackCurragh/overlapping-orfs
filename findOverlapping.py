from Bio import SeqIO
from Bio.Seq import Seq
from intervaltree import Interval, IntervalTree


def max_consecutive(arr):
    '''
    find the longest consecutive sequence of 0s in an array and location of start and end of sequence

    Input:
        arr: array of 0s and 1s

    Output:
        max_consecutive: integer of the longest consecutive sequence of 0s in the array
    '''
    consecutives = []
    current = [0,0]
    for i in range(len(arr)):
        if arr[i] == 0:
            current[1] = i
        else:
            consecutives.append(current)
            current = [i - 1, i]
    
    longest = [0, [0,0]]
    for start, end in consecutives:
        if end - start > longest[0]:
            longest = [end-start, [start, end]]
    
    return longest


def find_longest_overlapping_open_reading_length(mrna):
    '''
    Take mRNA sequence and idenfity the longest overlapping open reading frame (ORF) in the sequence.
    ORF is defined as the longest region of the mRNA sequence that does not contain a stop codon.

    Input: 
        mrna: string of mRNA sequence
    Output:
        longest_orf_length: integer of the length of the longest ORF
    
    '''
    stop_codons = ["TAA", "TAG", "TGA"]
    longest_orf_length = 0

    mrna_seq = Seq(mrna)

    # initialise array to store where stop codons are found. 0 for none 1 for found
    stop_codon_frames = [[],[],[]]
    for i in range(len(mrna_seq)):
        if mrna_seq[i:i+3] in stop_codons:
            stop_codon_frames[i%3].append(1)
        else:
            stop_codon_frames[i%3].append(0)

    longest_orf = 0
    longest_orfs = []
    frame1 = max_consecutive(stop_codon_frames[0])
    frame2 = max_consecutive(stop_codon_frames[1])
    frame3 = max_consecutive(stop_codon_frames[2])

    longest_orf = [
        frame1[0] * 3, frame1[1][0] * 3, frame1[1][1] * 3,
        frame2[0] * 3, frame2[1][0] * 3, frame2[1][1] * 3,
        frame3[0] * 3, frame3[1][0] * 3, frame3[1][1] * 3
    ]
    return longest_orf


def find_overlapping_regions(coordinates):
    '''
    Find the overlapping regions in a set of coordinates

    Input:
        coordinates: list of tuples of coordinates
    Output:
        overlapping_regions: list of tuples of overlapping regions
    '''
    # Build the Interval Tree
    interval_tree = IntervalTree(Interval(start, end) for start, end in coordinates)
    # Check overlaps in a pairwise manner and find 3-way overlaps with length 3
    three_way_overlaps = []
    overlaps = []

    for i, interval in enumerate(interval_tree):
        for other_interval in interval_tree:
            if interval == other_interval:  # Skip comparing an interval with itself
                continue
            if interval.overlaps(other_interval):
                overlap_start = max(interval.begin, other_interval.begin)
                overlap_end = min(interval.end, other_interval.end)
                overlap_length = max(0, overlap_end - overlap_start)
                if list(interval)[:2] not in [i[1] for i in overlaps]:
                    overlaps.append((list(interval)[:2], list(other_interval)[:2], overlap_length))
        
    print(overlaps)

def find_longest_regions_from_fasta(file_path):
    results = []

    for record in SeqIO.parse(file_path, "fasta"):
        mrna_name = record.id
        mrna_sequence = record.seq
        longest_region_length = find_longest_overlapping_open_reading_length(mrna_sequence)
        coordinates = [
            (longest_region_length[1], longest_region_length[2]),
            (longest_region_length[4], longest_region_length[5]),
            (longest_region_length[7], longest_region_length[8])
        ]
        find_overlapping_regions(coordinates)
        
        results.append((mrna_name, str(longest_region_length).strip("[]").replace(", ", "\t")))
        if len(results) > 2:
            return results

def write_results_to_tsv(results, output_file):
    with open(output_file, "w") as f:
        f.write("Sequence_Name\tFrame 0 Max\tFrame 0 Start\tFrame 0 End\tFrame 1 Max\tFrame 1 Start\tFrame 1 End\tFrame 2 Max\tFrame 2 Start\tFrame 2 End\n")
        for mrna_name, longest_length in results:
            f.write(f"{mrna_name}\t{longest_length}\n")


if __name__ == "__main__":
    # Example usage:
    fasta_file_path = "data/gencode.v44.transcripts.fa"
    output_tsv_file = "output_file.tsv"
    results = find_longest_regions_from_fasta(fasta_file_path)
    print(results)
    write_results_to_tsv(results, output_tsv_file)
