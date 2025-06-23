import os
import sys
import re
from typing import List, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import unicodedata

# Custom exceptions
class InvalidSequenceTypeError(Exception):
    """Raised when a non-DNA sequence (e.g., protein) is detected."""
    pass

class FileFormatError(Exception):
    """Raised when the file is not a valid FASTA file."""
    pass

class SequenceFormatError(Exception):
    """Raised when a sequence contains invalid DNA characters."""
    pass

class MissingSequenceError(Exception):
    """Raised when a header has no associated sequence."""
    pass

def validate_dna_sequence(sequence: str, seq_id: str) -> None:
    """Validate that the sequence is DNA and not protein or invalid."""
    valid_dna = set('ATGCN')
    sequence_upper = sequence.upper()
    
    # Check for invalid characters
    invalid_chars = set(sequence_upper) - valid_dna
    if invalid_chars:
        # Check if invalid characters are protein amino acids
        amino_acids = set('DEFHIKLMNPQRSVWY')
        if invalid_chars & amino_acids:
            raise InvalidSequenceTypeError(
                f"Sequence '{seq_id}' contains protein amino acids: {invalid_chars & amino_acids}."
            )
        raise SequenceFormatError(
            f"Sequence '{seq_id}' contains invalid DNA characters: {invalid_chars}"
        )

def clean_header(header: str, index: int) -> str:
    """Clean and standardize FASTA header."""
    # Remove leading '>' and strip whitespace
    header = header.lstrip('>').strip()
    if not header:
        header = f"sequence_{index}"
    # Remove control characters and normalize unicode
    header = ''.join(c for c in unicodedata.normalize('NFKD', header) if not unicodedata.category(c).startswith('C'))
    # Ensure header is not empty after cleaning
    if not header:
        header = f"sequence_{index}"
    return header

def format_fasta_sequence(sequence: str, line_length: int = 80) -> str:
    """Format sequence into lines of specified length."""
    return '\n'.join(sequence[i:i+line_length] for i in range(0, len(sequence), line_length))

def process_fasta_file(input_path: str, output_path: str) -> Tuple[int, List[str]]:
    """
    Read, validate, and clean FASTA file, producing an error-free output.
    Returns tuple of (number of sequences processed, list of warnings).
    """
    warnings = []
    valid_records = []
    seen_headers = set()
    
    try:
        # Check if input file exists and is readable
        if not os.path.exists(input_path):
            raise FileNotFoundError(f"Input file '{input_path}' does not exist.")
        if not os.access(input_path, os.R_OK):
            raise PermissionError(f"No read permission for file '{input_path}'.")
        if os.path.getsize(input_path) == 0:
            raise ValueError(f"Input file '{input_path}' is empty.")
        
        # Attempt to parse FASTA file
        try:
            records = list(SeqIO.parse(input_path, "fasta"))
        except ValueError as e:
            raise FileFormatError(f"Invalid FASTA format in '{input_path}': {str(e)}")
        
        if not records:
            raise ValueError(f"No sequences found in '{input_path}'.")
        
        for idx, record in enumerate(records, 1):
            seq_id = record.id or f"sequence_{idx}"
            sequence = str(record.seq).strip()
            
            try:
                # Check for missing sequence
                if not sequence:
                    warnings.append(f"Skipping header '{seq_id}': No sequence data.")
                    continue
                
                # Validate DNA sequence
                validate_dna_sequence(sequence, seq_id)
                
                # Clean and validate header
                cleaned_header = clean_header(record.description, idx)
                if cleaned_header in seen_headers:
                    warnings.append(f"Duplicate header '{cleaned_header}' detected; keeping both.")
                seen_headers.add(cleaned_header)
                
                # Create new sequence record
                valid_records.append(
                    SeqRecord(
                        Seq(sequence.upper()),
                        id=cleaned_header,
                        description=""
                    )
                )
            except (InvalidSequenceTypeError, SequenceFormatError, MissingSequenceError) as e:
                warnings.append(f"Skipping sequence '{seq_id}': {str(e)}")
                continue
        
        # Write cleaned FASTA file
        try:
            with open(output_path, 'w', encoding='utf-8') as output_handle:
                SeqIO.write(valid_records, output_handle, "fasta")
        except PermissionError:
            raise PermissionError(f"No write permission for output file '{output_path}'.")
        except IOError as e:
            raise IOError(f"Failed to write output file '{output_path}': {str(e)}")
        
        return len(valid_records), warnings
    
    except UnicodeDecodeError as e:
        raise UnicodeDecodeError(f"File '{input_path}' has invalid encoding: {str(e)}", e.object, e.start, e.end, e.reason)
    except MemoryError:
        raise MemoryError(f"Insufficient memory to process '{input_path}'.")
    except Exception as e:
        raise RuntimeError(f"Unexpected error processing '{input_path}': {str(e)}")

def main():
    """Command-line interface for the script."""
    parser = argparse.ArgumentParser(
        description="Clean and validate a FASTA file containing DNA gene sequences."
    )
    parser.add_argument(
        "input_file", help="Path to input FASTA file."
    )
    parser.add_argument(
        "output_file", help="Path to output cleaned FASTA file."
    )
    
    try:
        args = parser.parse_args()
    except SystemExit:
        raise ValueError("Invalid command-line arguments provided.")
    
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