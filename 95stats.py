import os
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter, defaultdict
import math
import re
from clean import process_fasta_file

def reverse_complement(s):
    return s.translate(str.maketrans('ATGC', 'TACG'))[::-1]

def calculate_tm_oligo(seq):
    """Calculate melting temperature for 20bp oligo using Wallace rule."""
    if len(seq) < 20:
        return "not applicable"
    sub = seq[:20]
    a = sub.count('A')
    t = sub.count('T')
    g = sub.count('G')
    c = sub.count('C')
    return f"{2 * (a + t) + 4 * (g + c)} °C"

def calculate_stats(seq_id, sequence):
    """Compute 95 statistics from a DNA sequence."""
    stats = {}
    seq = sequence.upper()
    length = len(seq)

    count = Counter(seq)
    A = count.get('A', 0)
    T = count.get('T', 0)
    G = count.get('G', 0)
    C = count.get('C', 0)
    N = count.get('N', 0)

    stats['1. Sequence Length'] = f"{length} bp"
    stats['2. A Count'] = A
    stats['3. T Count'] = T
    stats['4. G Count'] = G
    stats['5. C Count'] = C
    stats['6. N Count'] = N
    stats['7. A%'] = round(A / length * 100, 2) if length else 0
    stats['8. T%'] = round(T / length * 100, 2) if length else 0
    stats['9. G%'] = round(G / length * 100, 2) if length else 0
    stats['10. C%'] = round(C / length * 100, 2) if length else 0
    stats['11. N%'] = round(N / length * 100, 2) if length else 0
    stats['12. GC Content (%)'] = round((G + C) / length * 100, 2) if length else 0
    stats['13. AT Content (%)'] = round((A + T) / length * 100, 2) if length else 0
    stats['14. GC Skew'] = round((G - C) / (G + C), 3) if (G + C) else 'undefined (ratio)'
    stats['15. AT Skew'] = round((A - T) / (A + T), 3) if (A + T) else 'undefined (ratio)'

    dinuc_freq = Counter(seq[i:i+2] for i in range(len(seq) - 1))
    trinuc_freq = Counter(seq[i:i+3] for i in range(len(seq) - 2))
    for pair in sorted(dinuc_freq):
        stats[f"Dinucleotide {pair} Frequency"] = dinuc_freq[pair]
    for triplet in sorted(trinuc_freq):
        stats[f"Trinucleotide {triplet} Frequency"] = trinuc_freq[triplet]

    codons = [seq[i:i+3] for i in range(0, len(seq) - 2, 3)]
    codon_usage = Counter(codons)
    for codon in sorted(codon_usage):
        stats[f"Codon {codon} Count"] = codon_usage[codon]

    start_codon = 'ATG'
    stop_codons = {'TAA', 'TAG', 'TGA'}
    orfs = []
    i = 0
    while i < len(seq) - 2:
        codon = seq[i:i+3]
        if codon == start_codon:
            for j in range(i + 3, len(seq) - 2, 3):
                stop = seq[j:j+3]
                if stop in stop_codons:
                    orfs.append(seq[i:j+3])
                    break
        i += 3
    orf_lengths = [len(orf) for orf in orfs]
    stats['Start Codon Count'] = len(orfs)
    stats['Stop Codon Count'] = sum(seq.count(sc) for sc in stop_codons)
    stats['ORF Count'] = len(orfs)
    stats['Longest ORF Length'] = max(orf_lengths) if orf_lengths else 0
    stats['Shortest ORF Length'] = min(orf_lengths) if orf_lengths else 0

    motifs = ['TATAAA', 'AATAAA', 'GGGCGG', 'CCAAT']
    for motif in motifs:
        stats[f"Motif '{motif}' Count"] = seq.count(motif)

    cpg_islands = [m.start() for m in re.finditer(r'CG', seq)]
    stats['CpG Count'] = len(cpg_islands)

    homopolymers = re.findall(r'(A{6,}|T{6,}|G{6,}|C{6,})', seq)
    stats['Homopolymer Count (>6 bases)'] = len(homopolymers)

    reverse = seq[::-1]
    revcomp = str(Seq(seq).reverse_complement())
    stats['Reversed Sequence (start)'] = reverse[:20] + '...'
    stats['Reverse Complement (start)'] = revcomp[:20] + '...'

    hairpin_count = 0
    for stem_len in range(3, 5):
        for loop_len in range(3, 7):
            for i in range(len(seq) - 2*stem_len - loop_len):
                stem1 = seq[i:i+stem_len]
                loop = seq[i+stem_len:i+stem_len+loop_len]
                stem2 = seq[i+stem_len+loop_len:i+2*stem_len+loop_len]
                if reverse_complement(stem1) == stem2:
                    hairpin_count += 1
    stats['Hairpin Candidate Count'] = hairpin_count

    pattern_density = sum(len(set(seq[i:i+6] for i in range(len(seq) - 5))) for _ in range(1)) / length
    stats['Pattern Diversity Score'] = round(pattern_density, 4)

    prob = [n / length for n in [A, T, G, C] if n > 0]
    entropy = -sum(p * math.log2(p) for p in prob)
    stats['Shannon Entropy'] = round(entropy, 4)

    clamp = seq[-5:]
    gc_clamp_score = clamp.count('G') + clamp.count('C')
    stats['GC Clamp Score'] = f"{gc_clamp_score}/5 bases"

    stats['Melting Temperature (approx)'] = calculate_tm_oligo(seq)

    return stats

def process_sequences(input_file):
    clean_path = "__temp_cleaned__.fasta"
    _, warnings = process_fasta_file(input_file, clean_path)
    records = list(SeqIO.parse(clean_path, "fasta"))

    for record in records:
        stats = calculate_stats(record.id, str(record.seq))
        print(f"\n=== Stats for {record.id} ===")
        for key, value in stats.items():
            print(f"{key}: {value}")

    os.remove(clean_path)

def main():
    parser = argparse.ArgumentParser(description="Calculate 95 statistics from a DNA FASTA file.")
    parser.add_argument("input_file", help="Path to the input FASTA file.")
    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        print(f"File '{args.input_file}' not found.", file=sys.stderr)
        sys.exit(1)

    process_sequences(args.input_file)

if __name__ == "__main__":
    main()
