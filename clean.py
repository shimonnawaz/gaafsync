import os
import sys
import unicodedata
from typing import List, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
from detect_mistakes_in_fasta import validate_dna_sequence, clean_header, format_fasta_sequence, InvalidSequenceTypeError, SequenceFormatError, MissingSequenceError, FileFormatError

import logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Enable this to see debug messages
# logger.setLevel(logging.DEBUG)

def correct_sequence(sequence: str, seq_id: str) -> str:
    """Correct invalid DNA characters by replacing with 'N' and removing whitespace."""
    sequence = ''.join(sequence.split())
    sequence = ''.join(c if c.upper() in 'ATGCN' else 'N' for c in sequence)
    return sequence

def clean_and_unique_header(header: str, index: int, seen: set) -> str:
    """Clean header and ensure uniqueness."""
    base = clean_header(header, index)
    if base not in seen:
        return base
    count = 1
    new_header = f"{base}_{count}"
    while new_header in seen:
        count += 1
        new_header = f"{base}_{count}"
    return new_header

def process_fasta_file(input_path: str, output_path: str) -> Tuple[int, List[str]]:
    warnings = []
    summary = {"Skipped": 0, "Corrected": 0, "Duplicates": 0}
    valid_records = []
    seen_headers = set()

    try:
        if not os.path.exists(input_path):
            raise FileNotFoundError(f"Input file '{input_path}' does not exist.")
        if not os.access(input_path, os.R_OK):
            raise PermissionError(f"No read permission for file '{input_path}'.")
        if os.path.getsize(input_path) == 0:
            raise ValueError(f"Input file '{input_path}' is empty.")

        try:
            records = list(SeqIO.parse(input_path, "fasta"))
        except ValueError as e:
            raise FileFormatError(f"Invalid FASTA format in '{input_path}': {str(e)}")

        if not records:
            raise ValueError(f"No sequences found in '{input_path}'.")

        for idx, record in enumerate(records, 1):
            seq_id = record.id or f"sequence_{idx}"
            sequence = str(record.seq).strip()
            logger.debug(f"Processing sequence '{seq_id}' with data: {sequence}")

            if not sequence:
                warnings.append(f"Skipping header '{seq_id}': No sequence data.")
                summary["Skipped"] += 1
                continue

            cleaned_sequence = correct_sequence(sequence, seq_id)
            cleaned_header = clean_and_unique_header(record.description, idx, seen_headers)
            if cleaned_header != clean_header(record.description, idx):
                summary["Duplicates"] += 1

            seen_headers.add(cleaned_header)
            valid_records.append(
                SeqRecord(
                    Seq(cleaned_sequence.upper()),
                    id=cleaned_header,
                    description=""
                )
            )
            summary["Corrected"] += 1
            logger.debug(f"Added sequence '{seq_id}' as '{cleaned_header}'")

        try:
            with open(output_path, 'w', encoding='utf-8') as output_handle:
                for record in valid_records:
                    formatted_seq = format_fasta_sequence(str(record.seq))
                    output_handle.write(f">{record.id}\n{formatted_seq}\n")
        except PermissionError:
            raise PermissionError(f"No write permission for output file '{output_path}'.")
        except IOError as e:
            raise IOError(f"Failed to write output file '{output_path}': {str(e)}")

        # Append summary to warnings
        warnings.append(f"Summary: {summary['Skipped']} skipped, {summary['Corrected']} corrected, {summary['Duplicates']} renamed due to duplicates.")
        return len(valid_records), warnings

    except UnicodeDecodeError as e:
        raise UnicodeDecodeError(f"File '{input_path}' has invalid encoding: {str(e)}", e.object, e.start, e.end, e.reason)
    except MemoryError:
        raise MemoryError(f"Insufficient memory to process '{input_path}'.")
    except Exception as e:
        raise RuntimeError(f"Unexpected error processing '{input_path}': {str(e)}")

def main():
    parser = argparse.ArgumentParser(
        description="Clean and correct a FASTA file containing DNA gene sequences for data computation."
    )
    parser.add_argument("input_file", help="Path to input FASTA file.")
    parser.add_argument("output_file", help="Path to output cleaned FASTA file.")
    args = parser.parse_args()

    try:
        num_sequences, warnings = process_fasta_file(args.input_file, args.output_file)
        print(f"Successfully processed {num_sequences} sequences.")
        print(f"Output written to '{args.output_file}'.")
        if warnings:
            print("Warnings:")
            for warning in warnings:
                print(f"- {warning}")
    except (FileNotFoundError, PermissionError, ValueError, FileFormatError,
            UnicodeDecodeError, MemoryError, IOError, RuntimeError) as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
