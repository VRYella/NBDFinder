#!/usr/bin/env python3
"""
NBDFinder Batch Analysis Tool
============================

Command-line tool for batch processing of FASTA files.

Usage:
    python batch_analysis.py input.fasta output.csv [options]

Authors: Dr. Venkata Rajesh Yella
"""

import argparse
import pandas as pd
import sys
import time
from pathlib import Path
from Bio import SeqIO
from motifs import all_motifs

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="NBDFinder Batch Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python batch_analysis.py sequences.fasta results.csv
    python batch_analysis.py input.fa output.xlsx --format excel
    python batch_analysis.py seqs.fasta results.csv --motifs G4,Z-DNA
        """
    )
    
    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument("output", help="Output file (CSV or Excel)")
    parser.add_argument("--format", choices=["csv", "excel"], default="csv",
                       help="Output format (default: csv)")
    parser.add_argument("--motifs", help="Comma-separated list of motif types to analyze")
    parser.add_argument("--min-score", type=float, help="Minimum score threshold")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    
    return parser.parse_args()

def process_fasta_file(input_file, motif_filter=None, min_score=None, verbose=False):
    """Process FASTA file and return results"""
    results = []
    
    try:
        sequences = list(SeqIO.parse(input_file, "fasta"))
        total_sequences = len(sequences)
        
        if verbose:
            print(f"Processing {total_sequences} sequences...")
        
        for i, record in enumerate(sequences):
            if verbose and (i + 1) % 10 == 0:
                print(f"  Processed {i + 1}/{total_sequences} sequences")
            
            sequence_id = record.id
            sequence_str = str(record.seq).upper()
            
            # Analyze motifs
            motifs = all_motifs(sequence_str)
            
            # Apply filters
            filtered_motifs = motifs
            if motif_filter:
                filtered_motifs = [m for m in motifs if m.get('Class') in motif_filter]
            if min_score:
                filtered_motifs = [m for m in filtered_motifs 
                                 if m.get('Score', 0) >= min_score]
            
            # Add to results
            for motif in filtered_motifs:
                result = {
                    'Sequence_ID': sequence_id,
                    'Sequence_Length': len(sequence_str),
                    **motif  # Add all motif fields
                }
                results.append(result)
        
        if verbose:
            print(f"Found {len(results)} total motifs")
        
        return results
    
    except Exception as e:
        print(f"Error processing file: {e}")
        return []

def save_results(results, output_file, format_type="csv"):
    """Save results to file"""
    if not results:
        print("No results to save")
        return
    
    df = pd.DataFrame(results)
    
    try:
        if format_type == "excel":
            df.to_excel(output_file, index=False)
        else:
            df.to_csv(output_file, index=False)
        
        print(f"Results saved to {output_file}")
        print(f"Total motifs: {len(results)}")
        
        # Summary statistics
        motif_counts = df['Class'].value_counts()
        print("\nMotif summary:")
        for motif_type, count in motif_counts.items():
            print(f"  {motif_type}: {count}")
    
    except Exception as e:
        print(f"Error saving results: {e}")

def main():
    """Main function"""
    args = parse_arguments()
    
    # Validate input file
    if not Path(args.input).exists():
        print(f"Error: Input file '{args.input}' not found")
        sys.exit(1)
    
    # Parse motif filter
    motif_filter = None
    if args.motifs:
        motif_filter = [m.strip() for m in args.motifs.split(",")]
    
    # Process sequences
    start_time = time.time()
    results = process_fasta_file(
        args.input, 
        motif_filter=motif_filter,
        min_score=args.min_score,
        verbose=args.verbose
    )
    end_time = time.time()
    
    if args.verbose:
        print(f"Processing time: {end_time - start_time:.2f} seconds")
    
    # Save results
    save_results(results, args.output, args.format)

if __name__ == "__main__":
    main()
