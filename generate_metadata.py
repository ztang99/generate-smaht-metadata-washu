"""
generate_metadata.py

Orchestration script for generating comprehensive metadata sheets from sequencing CSV data.
Generates all required sheets for genomics data submission.

Author: Shadi Zaheri (refactored by Zitian Tang)
Date: 2025-05-23
Version: 2.0
"""

import os
import pandas as pd
import argparse
from pathlib import Path
from datetime import datetime

from sheet_generators import *
from utils import load_sequencing_data

def save_sheet(dataframe, tsv_path, xlsx_path, sheet_name):
    """Save dataframe to both TSV and Excel formats (if paths provided)."""
    if tsv_path:
        dataframe.to_csv(tsv_path, sep="\t", index=False)
        print(f"✅ {sheet_name} TSV saved to: {tsv_path}")
    if xlsx_path:
        dataframe.to_excel(xlsx_path, index=False)
        print(f"✅ {sheet_name} XLSX saved to: {xlsx_path}")
    if not tsv_path and not xlsx_path:
        print(f"✅ {sheet_name} sheet generated (in-memory only)")

def generate_all_metadata_sheets(args):
    """Generate all metadata sheets and save to specified outputs."""
    
    print("🔄 Loading sequencing data...")
    df = load_sequencing_data(args.input_csv)
    print(f"📊 Loaded {len(df)} sequencing records")
    
    print("🔄 Generating metadata sheets...")
    
    # ========================================================================
    # CORE METADATA SHEETS
    # ========================================================================
    
    # Donor Sheet - EMPTY (TPC handles this)
    print("📋 Generating Donor sheet...")
    donor_sheet = generate_donor_sheet(df, args.donor_prefix)
    save_sheet(donor_sheet, getattr(args, 'out_donor_tsv', None), getattr(args, 'out_donor_xlsx', None), "Donor")
    
    # Tissue Sheet - EMPTY (TPC handles this)
    print("📋 Generating Tissue sheet...")
    tissue_sheet = generate_tissue_sheet(df, args.donor_prefix)
    save_sheet(tissue_sheet, getattr(args, 'out_tissue_tsv', None), getattr(args, 'out_tissue_xlsx', None), "Tissue")
    
    # TissueSample Sheet
    print("📋 Generating TissueSample sheet...")
    tissue_sample_sheet = generate_tissue_sample_sheet(df, tissue_sheet, args.submitter_prefix, args.platform)
    save_sheet(tissue_sample_sheet, getattr(args, 'out_tissuesample_tsv', None), getattr(args, 'out_tissuesample_xlsx', None), "TissueSample")
    
    # Analyte Sheet
    print("📋 Generating Analyte sheet...")
    analyte_sheet = generate_analyte_sheet(df, tissue_sample_sheet, args.submitter_prefix, args.platform)
    save_sheet(analyte_sheet, getattr(args, 'out_analyte_tsv', None), getattr(args, 'out_analyte_xlsx', None), "Analyte")
    
    # AnalytePreparation Sheet
    print("📋 Generating AnalytePreparation sheet...")
    analyte_preparation_sheet = generate_analyte_preparation_sheet(df, args.submitter_prefix, args.platform)
    save_sheet(analyte_preparation_sheet, getattr(args, 'out_analytepreparation_tsv', None), getattr(args, 'out_analytepreparation_xlsx', None), "AnalytePreparation")
    
    # PreparationKit Sheet
    print("📋 Generating PreparationKit sheet...")
    preparation_kit_sheet = generate_preparation_kit_sheet(df, args.submitter_prefix, args.platform)
    save_sheet(preparation_kit_sheet, getattr(args, 'out_preparationkit_tsv', None), getattr(args, 'out_preparationkit_xlsx', None), "PreparationKit")
    
    # Treatment Sheet (empty)
    print("📋 Generating Treatment sheet...")
    treatment_sheet = generate_treatment_sheet()
    save_sheet(treatment_sheet, getattr(args, 'out_treatment_tsv', None), getattr(args, 'out_treatment_xlsx', None), "Treatment")
    
    # Library Sheet
    print("📋 Generating Library sheet...")
    library_sheet = generate_library_sheet(df, analyte_sheet, args.submitter_prefix, args.platform)
    save_sheet(library_sheet, getattr(args, 'out_library_tsv', None), getattr(args, 'out_library_xlsx', None), "Library")
    
    # LibraryPreparation Sheet
    print("📋 Generating LibraryPreparation sheet...")
    library_preparation_sheet = generate_library_preparation_sheet(args.submitter_prefix, args.platform)
    save_sheet(library_preparation_sheet, getattr(args, 'out_librarypreparation_tsv', None), getattr(args, 'out_librarypreparation_xlsx', None), "LibraryPreparation")
    
    # Sequencing Sheet
    print("📋 Generating Sequencing sheet...")
    sequencing_sheet = generate_sequencing_sheet(args.submitter_prefix, args.platform)
    save_sheet(sequencing_sheet, getattr(args, 'out_sequencing_tsv', None), getattr(args, 'out_sequencing_xlsx', None), "Sequencing")
    
    # ========================================================================
    # FILE-BASED SHEETS
    # ========================================================================
    
    # FileSet Sheet
    print("📋 Generating FileSet sheet...")
    fileset_sheet = generate_fileset_sheet(df, library_sheet, args.submitter_prefix, args.platform)
    save_sheet(fileset_sheet, getattr(args, 'out_fileset_tsv', None), getattr(args, 'out_fileset_xlsx', None), "FileSet")
    
    # UnalignedReads Sheet
    print("📋 Generating UnalignedReads sheet...")
    unaligned_reads_sheet = generate_unaligned_reads_sheet(df, fileset_sheet, args.submitter_prefix, args.platform)
    save_sheet(unaligned_reads_sheet, getattr(args, 'out_unalignedreads_tsv', None), getattr(args, 'out_unalignedreads_xlsx', None), "UnalignedReads")
    
    # AlignedReads Sheet (empty)
    print("📋 Generating AlignedReads sheet...")
    aligned_reads_sheet = generate_aligned_reads_sheet(df)
    save_sheet(aligned_reads_sheet, getattr(args, 'out_alignedreads_tsv', None), getattr(args, 'out_alignedreads_xlsx', None), "AlignedReads")
    
    # VariantCalls Sheet (empty)
    print("📋 Generating VariantCalls sheet...")
    variant_calls_sheet = generate_variant_calls_sheet()
    save_sheet(variant_calls_sheet, getattr(args, 'out_variantcalls_tsv', None), getattr(args, 'out_variantcalls_xlsx', None), "VariantCalls")
    
    # Software Sheet (empty)
    print("📋 Generating Software sheet...")
    software_sheet = generate_software_sheet()
    save_sheet(software_sheet, getattr(args, 'out_software_tsv', None), getattr(args, 'out_software_xlsx', None), "Software")

    # 250730 - NEW: SupplementaryFile Sheet (empty)
    print("📋 Generating SupplementaryFile sheet...")
    supplementary_file_sheet = generate_supplementary_file_sheet()
    save_sheet(supplementary_file_sheet, getattr(args, 'out_supplementaryfile_tsv', None), getattr(args, 'out_supplementaryfile_xlsx', None), "SupplementaryFile")
    
    
    # ========================================================================
    # COMBINED WORKBOOK
    # ========================================================================
    
    print("📊 Creating combined Excel workbook...")
    create_combined_workbook(args, {
        'donor': donor_sheet,
        'tissue': tissue_sheet,
        'tissue_sample': tissue_sample_sheet,
        'analyte': analyte_sheet,
        'analyte_preparation': analyte_preparation_sheet,
        'preparation_kit': preparation_kit_sheet,
        'treatment': treatment_sheet,
        'library': library_sheet,
        'library_preparation': library_preparation_sheet,
        'sequencing': sequencing_sheet,
        'fileset': fileset_sheet,
        'unaligned_reads': unaligned_reads_sheet,
        'aligned_reads': aligned_reads_sheet,
        'variant_calls': variant_calls_sheet,
        'supplementary_file': supplementary_file_sheet,
        'software': software_sheet
    })

def create_combined_workbook(args, sheets_dict):
    """Create combined Excel workbook with all metadata sheets."""
    with pd.ExcelWriter(args.out_combined_xlsx, engine="openpyxl") as writer:
        # Add overview guidelines sheet (placeholder)
        overview_df = pd.DataFrame()
        overview_df.to_excel(writer, sheet_name="(Overview Guidelines)", index=False)
        
        # Add all metadata sheets
        sheets_dict['donor'].to_excel(writer, sheet_name="Donor", index=False)
        sheets_dict['tissue'].to_excel(writer, sheet_name="Tissue", index=False)
        sheets_dict['tissue_sample'].to_excel(writer, sheet_name="TissueSample", index=False)
        sheets_dict['analyte'].to_excel(writer, sheet_name="Analyte", index=False)
        sheets_dict['analyte_preparation'].to_excel(writer, sheet_name="AnalytePreparation", index=False)
        sheets_dict['preparation_kit'].to_excel(writer, sheet_name="PreparationKit", index=False)
        sheets_dict['treatment'].to_excel(writer, sheet_name="Treatment", index=False)
        sheets_dict['library'].to_excel(writer, sheet_name="Library", index=False)
        sheets_dict['library_preparation'].to_excel(writer, sheet_name="LibraryPreparation", index=False)
        sheets_dict['sequencing'].to_excel(writer, sheet_name="Sequencing", index=False)
        sheets_dict['fileset'].to_excel(writer, sheet_name="FileSet", index=False)
        sheets_dict['unaligned_reads'].to_excel(writer, sheet_name="UnalignedReads", index=False)
        sheets_dict['aligned_reads'].to_excel(writer, sheet_name="AlignedReads", index=False)
        sheets_dict['variant_calls'].to_excel(writer, sheet_name="VariantCalls", index=False)
        sheets_dict['supplementary_file'].to_excel(writer, sheet_name="SupplementaryFile", index=False)  # 250730-NEW
        sheets_dict['software'].to_excel(writer, sheet_name="Software", index=False)
    
    print(f"✅ Combined metadata Excel workbook saved to: {args.out_combined_xlsx}")

def setup_argument_parser():
    """Setup command line argument parser."""
    parser = argparse.ArgumentParser(
        description="Generate comprehensive metadata sheets from sequencing CSV data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
            Examples:
            # Generate only combined Excel workbook
            python generate_metadata.py \\
                --input-csv sequencing_data.csv \\
                --donor-prefix NDRI \\
                --submitter-prefix WASHU \\
                --out-combined-xlsx output/metadata_workbook.xlsx
                
            # Generate combined workbook plus individual files
            python generate_metadata.py \\
                --input-csv sequencing_data.csv \\
                --donor-prefix NDRI \\
                --submitter-prefix WASHU \\
                --out-donor-tsv output/donor.tsv \\
                --out-donor-xlsx output/donor.xlsx \\
                [... other output paths ...] \\
                --out-combined-xlsx output/metadata_workbook.xlsx
        """
    )
    
    # Input files
    input_group = parser.add_argument_group('Input Files')

    input_group.add_argument(
        "--work-dir", required=True,
        help="Working directory containing input CSV and where output will be saved"
    )
    input_group.add_argument(
        "--in-csv", required=True,
        help="Filename of the input CSV file (located in work-dir)"
    )
    
    # Configuration
    config_group = parser.add_argument_group('Configuration')
    config_group.add_argument(
        "--donor-prefix", default="NDRI",
        help="Prefix for donor and tissue IDs (default: NDRI)"
    )
    config_group.add_argument(
        "--submitter-prefix", default="WASHU",
        help="Prefix for other submitted_id fields (default: WASHU)"
    )

    config_group.add_argument(
        "--platform", default="Illumina", choices=["Illumina", "RNAseq", "PacBio"],
        help="Sequencing platform (default: Illumina). Options: Illumina, RNAseq, PacBio"
    )
    
    # Output files - TSV (optional)
    tsv_group = parser.add_argument_group('TSV Output Files (Optional)')
    tsv_outputs = [
        "donor", "tissue", "tissuesample", "analyte", "analytepreparation", 
        "preparationkit", "treatment", "library", "librarypreparation", 
        "sequencing", "fileset", "unalignedreads", "alignedreads", 
        "variantcalls", "software"
    ]
    for output in tsv_outputs:
        tsv_group.add_argument(
            f"--out-{output}-tsv", required=False,
            help=f"{output.replace('_', ' ').title()} sheet TSV output path (optional)"
        )
    
    # Output files - XLSX (optional)
    xlsx_group = parser.add_argument_group('XLSX Output Files (Optional)')
    for output in tsv_outputs:
        xlsx_group.add_argument(
            f"--out-{output}-xlsx", required=False,
            help=f"{output.replace('_', ' ').title()} sheet XLSX output path (optional)"
        )
    
    # Combined output
    ## 2025/10/30: Now will be auto-generated based on work-dir, sample name, platform, and date
    # parser.add_argument(
    #     "--out-combined-xlsx", required=True,
    #     help="Combined Excel workbook containing all sheets"
    # )
    
    return parser

def validate_inputs(args):
    """Validate input files exist and are readable."""
    # Construct full input path
    input_csv_path = Path(args.work_dir) / args.in_csv
    if not input_csv_path.exists():
        raise FileNotFoundError(f"Input CSV file not found: {input_csv_path}")
    if not input_csv_path.is_file():
        raise ValueError(f"Path is not a file: {input_csv_path}")
    
    # Validate work directory exists
    if not Path(args.work_dir).exists():
        raise FileNotFoundError(f"Work directory not found: {args.work_dir}")

def setup_file_paths(args):
    """Setup input and output file paths based on arguments."""
    # Construct input CSV path
    args.input_csv = str(Path(args.work_dir) / args.in_csv)
    
    # Get info to generate output filename
    sample_name = Path(args.work_dir).name
    current_date = datetime.now().strftime("%Y%m%d")
    
    output_filename = f"{sample_name}_{args.platform}_{current_date}.xlsx"
    
    # Set the combined output path
    args.out_combined_xlsx = str(Path(args.work_dir) / output_filename)
    
    print(f"📁 Work directory: {args.work_dir}")
    print(f"📄 Input CSV: {args.input_csv}")
    print(f"📊 Output file: {args.out_combined_xlsx}")

def main():
    """Main entry point for metadata generation."""
    parser = setup_argument_parser()
    args = parser.parse_args()
    
    try:
        # Validate inputs
        print("🔍 Validating input files...")
        validate_inputs(args)
        
        # Setup file paths
        print("🔧 Setting up file paths...")
        setup_file_paths(args)

        # Create output directories if needed
        output_paths = [args.out_combined_xlsx]
        
        # Add individual output paths if provided
        for attr_name in dir(args):
            if attr_name.startswith('out_') and attr_name.endswith(('_tsv', '_xlsx')):
                output_path = getattr(args, attr_name)
                if output_path:
                    output_paths.append(output_path)
        
        for output_path in output_paths:
            Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        
        # Generate all metadata sheets
        generate_all_metadata_sheets(args)
        
        print("🎉 Metadata generation completed successfully!")
        
    except Exception as e:
        print(f"❌ Error during metadata generation: {e}")
        raise

if __name__ == "__main__":
    main()