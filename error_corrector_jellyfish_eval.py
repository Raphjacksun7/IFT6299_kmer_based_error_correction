import os
import subprocess
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict
import numpy as np
import csv


def reverse_complement(kmer):
    """
    Generate the reverse complement of a DNA k-mer.
    
    Args:
        kmer (str): DNA sequence to be reversed and complemented.
    
    Returns:
        str: The reverse complement of the k-mer.
    """
    complement = str.maketrans('ACGT', 'TGCA')
    return kmer.translate(complement)[::-1]



def initialize_weights(k):
    """
    Initializes weights for k-mer positions and nucleotide bases (A, C, G, T) to 1.0.
    
    Args:
        k (int): Length of the k-mer.
    
    Returns:
        numpy.ndarray: Array with shape (k, 4) filled with 1.0.
    """
    # Create an array of ones with shape (k, 4) for unbiased initial weights
    return np.ones((k, 4))


def adjust_weights(weights, position, base_index, success):
    """
    Adjusts weights based on correction success, increasing for success and decreasing for failure, then normalizes.
    
    Args:
        weights (numpy.ndarray): Current weight matrix (k, 4).
        position (int): Position within the k-mer for adjustment.
        base_index (int): Index of the base (0: 'A', 1: 'C', 2: 'G', 3: 'T').
        success (bool): True if the correction succeeded, False otherwise.
    
    Modifies:
        weights (numpy.ndarray): Adjusted in-place to reflect correction success or failure.
    """
    # Set adjustment factor based on correction outcome
    adjustment_factor = 1.1 if success else 0.95
    # Apply adjustment to the specific weight and normalize the row
    weights[position, base_index] *= adjustment_factor
    weights[position] /= weights[position].sum()


def select_base(position, weights, current_base):
    """
    Selects a base for correction using a weighted approach with an epsilon-greedy strategy to avoid local minima.

    Args:
        position (int): Position in the k-mer for potential correction.
        weights (numpy.ndarray): Matrix of weights for each base at each position.
        current_base (str): Current base at the specified position.

    Returns:
        str: New base selected for correction. Returns current base if it's non-standard.
    """
    bases = "ACGT"  # Define the standard DNA bases.
    
    # Return the current base unchanged if it is not one of the standard bases (e.g., 'N').
    if current_base not in bases:
        return current_base

    # Identify the index of the current base and prepare alternative choices.
    idx = bases.index(current_base)
    choices = np.delete(np.arange(4), idx)  # Exclude the current base from choices.

    # Epsilon-greedy strategy:
    # With 5% probability, choose a random base to encourage exploration.
    if np.random.random() < 0.05:  
        return bases[np.random.choice(choices)]
    
    # With 95% probability, choose the base with the highest weight, promoting exploitation.
    return bases[choices[np.argmax(weights[position, choices])]]



def correct_read(read, kmer_counts, k, weights, quality_scores, i0):
    """
    Corrects sequencing errors in a given DNA read by examining each k-mer and potentially adjusting one base per k-mer
    based on k-mer frequency and base quality scores using a weighted selection strategy.

    Args:
        read (str): The DNA sequence to be corrected.
        kmer_counts (dict): A dictionary containing the counts of each k-mer and its reverse complement.
        k (int): The length of the k-mers to consider within the read.
        weights (numpy.ndarray): The weight matrix used for deciding the most probable base at each position.
        quality_scores (list): A list of quality scores corresponding to each base in the read.
        i0 (int): The threshold below which a k-mer is considered rare.

    Returns:
        tuple: A tuple containing the corrected DNA sequence as a string and the number of corrections made.
    """
    
    corrected_read = list(read)  # Convert the immutable string to a mutable list.
    num_corrections = 0  # Initialize count of successful corrections.
    
    for i in range(len(read) - k + 1):
        kmer = read[i:i+k]
        rc_kmer = reverse_complement(kmer)
        combined_count = kmer_counts.get(kmer, 0) + kmer_counts.get(rc_kmer, 0)
        
        # Check if k-mer is rare.
        if combined_count < i0:
            # Find index of the base with the lowest quality score within the k-mer
            min_q_index = min(range(k), key=lambda x: quality_scores[i + x])
            min_q_base = kmer[min_q_index]
            
            # Skip correction if the base is non-standard
            if min_q_base not in "ACGT":
                continue

            best_base = select_base(min_q_index, weights, min_q_base) # Select the best replacement base.

            if best_base != min_q_base:
                new_kmer = kmer[:min_q_index] + best_base + kmer[min_q_index+1:]
                new_rc_kmer = reverse_complement(new_kmer)
                new_combined_count = kmer_counts.get(new_kmer, 0) + kmer_counts.get(new_rc_kmer, 0)
                
                # Apply correction if the new k-mer is more frequent.
                if new_combined_count > i0:
                    corrected_read[i + min_q_index] = best_base
                    adjust_weights(weights, min_q_index, "ACGT".index(best_base), True) # Positively adjust weights.
                    num_corrections += 1
                    break  # Correct only one base per k-mer
                else:
                    adjust_weights(weights, min_q_index, "ACGT".index(min_q_base), False) # Negatively adjust weights.
    return "".join(corrected_read), num_corrections

def compute_frequency_spectrum(kmer_counts):
    """
    Computes the frequency spectrum of k-mer counts.
    Args:
        kmer_counts (dict): K-mers Dictionary.
    Returns:
        defaultdict: A dictionary mapping each k-mer count.
    """
    
    frequency_spectrum = defaultdict(int)
    for count in kmer_counts.values():
        frequency_spectrum[count] += 1
    return frequency_spectrum

def find_i0_and_imax(frequency_spectrum):
    """
    Finds the first minimum (i0) and the first maximum (imax) in the frequency spectrum after i0.
    Args:
        frequency_spectrum (dict): Frequency spectrum of k-mer counts.
    Returns:
        tuple: (i0, imax)
    """
    
    sorted_items = sorted(frequency_spectrum.items())
    i0, imax = None, None
    for i in range(1, len(sorted_items)-1):
        # Check for the first minimum.
        if not i0 and sorted_items[i-1][1] > sorted_items[i][1] < sorted_items[i+1][1]:
            i0 = sorted_items[i][0]
        # Check for the first maximum after the first minimum.
        elif i0 and sorted_items[i-1][1] < sorted_items[i][1] > sorted_items[i+1][1]:
            imax = sorted_items[i][0]
            break  # Stop searching after finding the first maximum.
    return i0, imax

def run_jellyfish(input_file, k, output_file):
    """
    Runs Jellyfish to count k-mers in a given FASTQ file.
    
    Args:
        input_file (str): Path to the input FASTQ file.
        k (int): The k-mer length to use.
        output_file (str): Path where the Jellyfish output should be saved.
    
    Raises:
        RuntimeError: If Jellyfish fails to run correctly.
    """
    try:
        subprocess.run([
            "jellyfish", "count", "-m", str(k), "-s", "500M", "-t", "8", "-C",
            "-o", output_file, input_file
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running Jellyfish count: {e}")
        raise RuntimeError("Failed to count k-mers with Jellyfish.")

    
def dump_jellyfish_kmers_to_file(jellyfish_file, output_dump):
    """
    Dumps the k-mer counts from Jellyfish output to a readable text file.
    
    Args:
        jellyfish_file (str): Path to the Jellyfish output file.
        output_dump (str): Path where the dump should be saved.
    
    Raises:
        RuntimeError: If the dump operation fails.
    """
    try:
        subprocess.run([
            "jellyfish", "dump", "-c", jellyfish_file, "-o", output_dump
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error dumping k-mers from Jellyfish: {e}")
        raise RuntimeError("Failed to dump k-mers with Jellyfish.")

def load_kmer_counts(dump_file):
    """
    Loads k-mer counts from a Jellyfish dump file into a dictionary.

    Args:
        dump_file (str): Path to the file containing dumped k-mer counts.

    Returns:
        defaultdict: K-mers dictionary.
    """
    kmer_counts = defaultdict(int)
    with open(dump_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            kmer = parts[0]
            count = int(parts[1])
            kmer_counts[kmer] = count
    return kmer_counts



def calculate_metrics(kmer_counts, i0):
    """
    Calculates metrics based on k-mer counts: number of rare k-mers and total k-mers.

    Args:
        kmer_counts (dict): Dictionary of k-mer counts.
        i0 (int): Threshold for determining if a k-mer is considered rare.

    Returns:
        tuple: A tuple containing the number of rare k-mers and the total number of k-mers.
    """
    rare_kmers = sum(1 for count in kmer_counts.values() if count <= i0)
    total_kmers = sum(kmer_counts.values())
    return rare_kmers, total_kmers # Return the calculated metrics as a tuple.

def write_metrics_to_csv(before_metrics, after_metrics, output_file, num_corrections):
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        writer.writerow(["Metric", "Before Correction", "After Correction", "Change %"])
        metrics = ["Rare K-mers", "Total K-mer Occurrences"]
        for metric, before, after in zip(metrics, before_metrics, after_metrics):
            change_percent = ((after - before) / before * 100) if before else 0
            writer.writerow([metric, before, after, f"{change_percent:.2f}%"])
            
        writer.writerow(["Corrections Made", "- ", num_corrections, "--"])
            
            


def perform_evaluation(output_count_file, k, output_directory, kmer_counts, weights, correction_ratio, total_corrections, i0):
    try:
        # File paths for the new Jellyfish results
        jellyfish_corrected_result_file = os.path.join(output_directory, "eval_kmer_counts.jf")
        jellyfish_corrected_dump_file = os.path.join(output_directory, "eval_kmer_counts_dump.txt")

        # Run Jellyfish on the corrected FASTQ file
        run_jellyfish(output_count_file, k, jellyfish_corrected_result_file)
        dump_jellyfish_kmers_to_file(jellyfish_corrected_result_file, jellyfish_corrected_dump_file)

        # Load corrected k-mer counts
        corrected_kmer_counts = load_kmer_counts(jellyfish_corrected_dump_file)

        # Calculate metrics
        original_rare_kmers, original_total_kmers = calculate_metrics(kmer_counts, i0)
        corrected_rare_kmers, corrected_total_kmers = calculate_metrics(corrected_kmer_counts, i0)
        
        # Compute final weights for each nucleotide
        final_weights = np.sum(weights, axis=0) / np.sum(weights)
        
        # Print the final weights
        print("Final Weights:")
        print(final_weights)

        # Generate log-like output
        eval_results_file = os.path.join(output_directory, "eval_results.txt")
        with open(eval_results_file, 'w') as logfile:
            logfile.write(f"Evaluation Results:\n")
            logfile.write(f"Metric: Total Corrections Made: {total_corrections:,}\n")
            logfile.write(f"Metric: Corrected vs. Uncorrected Read Ratio:  {correction_ratio:.2f}% \n")
            logfile.write(f"Metric: Rare K-mers | Before Correction: {original_rare_kmers:,} | After Correction: {corrected_rare_kmers:,} | Variation: {(corrected_rare_kmers - original_rare_kmers) / original_rare_kmers * 100:.4f}%\n")
            logfile.write(f"Final Weights:\n")
            for base, weight in zip("ACGT", final_weights):
                logfile.write(f"{base}: {weight:.4f}\n")

        print("Evaluation completed. Results saved to:", eval_results_file)
    except Exception as e:
        print(f"Error during evaluation: {e}")





def main(input_file, k, output_directory, evaluate=False, max_reads=None):
    try:
        base_filename = os.path.basename(input_file)
        output_count_file = os.path.join(output_directory, f"{os.path.splitext(base_filename)[0]}_corrected.fastq")
        
        os.makedirs(output_directory, exist_ok=True)
        jellyfish_result_file = os.path.join(output_directory, "kmer_counts.jf")
        jellyfish_dump_file = os.path.join(output_directory, "kmer_counts_dump.txt")
        
        run_jellyfish(input_file, k, jellyfish_result_file)
        dump_jellyfish_kmers_to_file(jellyfish_result_file, jellyfish_dump_file)
        kmer_counts = load_kmer_counts(jellyfish_dump_file)

        frequency_spectrum = compute_frequency_spectrum(kmer_counts)
        i0, _ = find_i0_and_imax(frequency_spectrum)

        weights = initialize_weights(k)
        reads_with_corrections = 0
        corrected_records = []
        total_corrections = 0
        read_count = 0

        for record in SeqIO.parse(input_file, "fastq"):
            if max_reads is not None and read_count >= max_reads:
                break
            corrected_seq, num_corrections = correct_read(str(record.seq), kmer_counts, k, weights, record.letter_annotations['phred_quality'], i0)
            corrected_records.append(SeqRecord(Seq(corrected_seq), id=record.id, description=record.description, letter_annotations={'phred_quality': record.letter_annotations['phred_quality']}))
            total_corrections +=num_corrections
            read_count += 1
            if num_corrections > 0:
                reads_with_corrections += 1

        with open(output_count_file, "w") as output_handle:
            SeqIO.write(corrected_records, output_handle, "fastq")

        print("Error correction completed and output saved.")

        if evaluate:
            correction_ratio = reads_with_corrections / read_count * 100 if read_count > 0 else 0
            perform_evaluation(output_count_file, k, output_directory, kmer_counts, weights, correction_ratio, total_corrections, i0)

    except Exception as e:
        print(f"An error occurred: {e}")
        
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="K-mer based error correction tool.")
    parser.add_argument("input_file", help="Input FASTQ file")
    parser.add_argument("-k", type=int, help="K-mer length", default=17)
    parser.add_argument("-d", "--directory", help="Output directory for results", default="output")
    parser.add_argument("--evaluate", action='store_true', help="Enable evaluation of correction efficiency")
    parser.add_argument("--max_reads", type=int, help="Maximum number of reads to process", default=None)
    args = parser.parse_args()

    main(args.input_file, args.k, args.directory, args.evaluate, args.max_reads)