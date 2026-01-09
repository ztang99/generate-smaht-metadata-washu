"""
strip_file_paths.py

Utility to strip full paths from file columns, keeping only filenames.
Processes 'Read1 File' and 'Read2 File' columns in CSV files.

Author: Shadi Zaheri (refactored by Zitian Tang)
Date: 2025-05-23
"""

import csv
import argparse
from pathlib import Path

def strip_file_paths(input_file, output_file, columns_to_strip=None):
    """
    Strip full paths from specified columns, keeping only filenames.
    
    Args:
        input_file: Path to input CSV file
        output_file: Path to output CSV file
        columns_to_strip: List of column names to process (default: ['Read1 File', 'Read2 File'])
    """
    if columns_to_strip is None:
        columns_to_strip = ['Read1 File', 'Read2 File']
    
    rows_processed = 0
    paths_stripped = 0
    
    with open(input_file, "r") as infile, open(output_file, "w", newline="") as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        
        writer.writeheader()
        
        for row in reader:
            rows_processed += 1
            
            # Process each specified column
            for column in columns_to_strip:
                if column in row and row[column]:
                    original_path = row[column]
                    filename = Path(original_path).name
                    
                    if original_path != filename:  # Only count if actually stripped
                        paths_stripped += 1
                        row[column] = filename
            
            writer.writerow(row)
    
    print(f"✅ Processed {rows_processed} rows")
    print(f"📁 Stripped {paths_stripped} file paths")
    print(f"💾 Output written to: {output_file}")

def main():
    """Main entry point for path stripping utility."""
    parser = argparse.ArgumentParser(
        description="Strip full paths from file columns, keeping only filenames",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Strip default columns (Read1 File, Read2 File)
  python strip_file_paths.py --input data.csv --output data_stripped.csv
  
  # Strip custom columns
  python strip_file_paths.py --input data.csv --output data_stripped.csv --columns "Read1 File" "Read2 File"
        """
    )
    
    parser.add_argument(
        "--input", required=True,
        help="Path to input CSV file"
    )
    parser.add_argument(
        "--output", required=True,
        help="Path to save output CSV with stripped paths"
    )
    parser.add_argument(
        "--columns", nargs="+", default=['Read1 File', 'Read2 File'],
        help="Column names to strip paths from (default: 'Read1 File' 'Read2 File')"
    )
    
    args = parser.parse_args()
    
    try:
        # Validate input file exists
        if not Path(args.input).exists():
            raise FileNotFoundError(f"Input file not found: {args.input}")
        
        # Create output directory if needed
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        
        # Strip paths
        strip_file_paths(args.input, args.output, args.columns)
        
    except Exception as e:
        print(f"❌ Error: {e}")
        raise

if __name__ == "__main__":
    main()