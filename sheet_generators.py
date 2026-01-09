"""
sheet_generators.py

Individual sheet generation functions for metadata workbook.
Each function creates a specific sheet type (Donor, Tissue, etc.)
"""

import pandas as pd
from utils import *

# ============================================================================
# CORE METADATA SHEETS
# ============================================================================

def generate_donor_sheet(df, prefix="NDRI"):
    """Generate donor metadata sheet from sequencing data."""
    return pd.DataFrame(columns=[
        'submitted_id', 'age', 'external_id', 'sex', 'tpc_submitted', 
        'eligibility', 'hardy_scale', 'protected_donor'  # 250730-added protected_donor, reordered
    ])


def generate_tissue_sheet(df, prefix="NDRI"):
    """Generate tissue metadata sheet from sequencing data."""
    return pd.DataFrame(columns=[
        'submitted_id', 'external_id', 'sample_count', 'preservation_type', 
        'preservation_medium', 'ischemic_time', 'anatomical_location', 
        'pathology_notes', 'ph', 'prosector_notes', 'size', 'size_unit', 
        'volume', 'weight', 'donor', 'uberon_id'  # 250730 - uberon_id moved to end
    ])


def generate_tissue_sample_sheet(df, tissue_sheet, prefix="WASHU", platform="Illumina"):
    """Generate tissue sample metadata sheet from sequencing data."""
    sample_df = get_unique_samples(df)
    tissue_groups = sample_df.groupby(['donor_id', 'tissue_code']).first().reset_index()
    
    processing_notes_map = {
        'Illumina': 'DNeasy',
        'RNAseq': 'RNeasy', 
        'PacBio': 'MagAttract'
    }
    
    def get_processing_notes(tissue_name, platform):
        if tissue_name == 'BLOOD' and platform == 'Illumina':
            return 'Promega Blood Kit'
        if platform not in processing_notes_map:
            raise ValueError(f"Unsupported platform: {platform}. Supported platforms: {list(processing_notes_map.keys())}")
        return processing_notes_map[platform]
    
    sample_records = []
    for _, row in tissue_groups.iterrows():
        processing_notes = get_processing_notes(row['tissue_name'], platform)
        # Generate sample_sources using direct mapping
        sample_sources = generate_sample_sources_id(row['donor_id'], row['tissue_code'], donor_prefix="NDRI")
        
        sample_records.append({
            'submitted_id': format_tissue_sample_id(prefix, row['sample_name']),
            'category': determine_category(row['tissue_name']),
            'external_id': row['sample_name'],
            'preservation_type': '',
            'preservation_medium': '',
            'description': '',
            'core_size': '',
            'processing_date': '',
            'processing_notes': processing_notes,
            'weight': '',
            'sample_sources': sample_sources,
            'parent_samples': ''
        })
    
    return pd.DataFrame(sample_records)


def generate_analyte_sheet(df, tissue_sample_sheet, prefix="WASHU", platform="Illumina"):
    """Generate analyte metadata sheet from sequencing data."""
    sample_df = get_unique_samples(df)
    tissue_groups = sample_df.groupby(['donor_id', 'tissue_code']).first().reset_index()
    ts_map = dict(zip(tissue_sample_sheet['external_id'], tissue_sample_sheet['submitted_id']))
    
    if platform == 'RNAseq':
        molecule = 'RNA'
        molecule_detail = 'Total RNA'
        suffix = '_GRNA'
    else:  # Illumina or PacBio
        molecule = 'DNA'
        molecule_detail = 'Total DNA'
        suffix = '_GDNA'
    
    analyte_records = []

    for _, row in tissue_groups.iterrows():
        kit_name = get_preparation_kit_for_tissue(row['tissue_name'], platform)
        
        tissue_sample_id = ts_map.get(row['sample_name'], '')
        if tissue_sample_id:
            analyte_submitted_id = tissue_sample_id.replace('TISSUE-SAMPLE', 'ANALYTE') + suffix
        else:
            analyte_submitted_id = format_analyte_id(prefix, row['sample_name'], suffix)
        
        # Add all required columns
        analyte_records.append({
            'submitted_id': analyte_submitted_id,
            'molecule': molecule,
            'molecule_detail': molecule_detail,
            'external_id': row['sample_name'],
            'description': '',
            'a260_a280_ratio': '',
            'average_fragment_size': '',
            'concentration': '',
            'concentration_unit': '', # 20250804 no concentration unit we don't need to fill this out
            'dna_integrity_number': '',
            'dna_integrity_number_instrument': '',
            'dna_quality_number': '',
            'dna_quality_number_instrument': '',
            'dna_quality_size_threshold': '',
            'genomic_quality_number': '',
            'genomic_quality_number_instrument': '',
            'genomic_quality_size_threshold': '',
            'quantitation_method': '',
            'ribosomal_rna_ratio': '',
            'rna_integrity_number': '',
            'rna_integrity_number_instrument': 'Agilent 4200 TapeStation' if platform == 'RNAseq' else '',
            'sample_quantity': '',
            'sample_quantity_unit': '',
            'volume': '',
            'volume_unit': '',
            'total_yield': '',
            'yield_per_unit': '',
            'yield_unit': '',
            'samples': ts_map.get(row['sample_name'], ''),
            'analyte_preparation': format_analyte_preparation_id(prefix, kit_name)
        })
    
    return pd.DataFrame(analyte_records)


def generate_analyte_preparation_sheet(df, prefix="WASHU", platform="Illumina"):
    """Generate analyte preparation metadata sheet."""
    sample_df = get_unique_samples(df)
    tissue_types = sample_df['tissue_name'].unique()

    prep_details = {}
    
    if platform == 'RNAseq':
        prep_details['RNEASY_KIT_PLUS'] = {
            'description': 'RNeasy Plus Universal Mini Kit',
            'cell_lysis_method': '',
            'extraction_method': ''
        }
    else:  # Illumina or PacBio
        if 'BLOOD' in tissue_types and platform == 'Illumina':
            prep_details['PROMEGA_RSC_BLOOD_DNA_KIT'] = {
                'description': '',
                'cell_lysis_method': '',
                'extraction_method': ''
            }
        
        if platform == 'Illumina':
            prep_details['DNEASY_KIT'] = {
                'description': 'DNeasy Blood and Tissue Kit',
                'cell_lysis_method': '',
                'extraction_method': ''
            }
        else:  # PacBio or mixed
            prep_details['MAGATTRACT_HMW_DNA_KIT'] = {
                'description': '',
                'cell_lysis_method': 'Enzymatic',
                'extraction_method': 'Magnetic Beads '
            }
    
    prep_records = []
    for kit_name, details in prep_details.items():
        prep_records.append({
            'submitted_id': format_analyte_preparation_id(prefix, kit_name),
            'description': details['description'],
            'cell_lysis_method': details['cell_lysis_method'],
            'cell_selection_method': '',  # 250730 - CHANGED: was cell_sorting_method
            'extraction_method': details['extraction_method'],
            'homogenization_method': '',
            'suspension_type': '',
            'preparation_kits': format_preparation_kit_id(prefix, kit_name),
            'treatments': ''
        })
    
    return pd.DataFrame(prep_records)


def generate_preparation_kit_sheet(df, prefix="WASHU", platform="Illumina"):
    """Generate preparation kit metadata sheet."""
    
    sample_df = get_unique_samples(df)
    tissue_types = sample_df['tissue_name'].unique()
    
    kit_details = {}
    
    if platform == 'RNAseq':
        kit_details['RNEASY_KIT_PLUS'] = {
            'title': 'RNeasy Kit Plus',
            'vendor': 'Qiagen',
            'version': '',
            'catalog_number': '73404'
        }
    else:  # Illumina or PacBio
        if 'BLOOD' in tissue_types and platform == 'Illumina':
            kit_details['PROMEGA_RSC_BLOOD_DNA_KIT'] = {
                'title': 'Promega RSC Blood DNA Kit',
                'vendor': 'Promega',
                'version': '',
                'catalog_number': 'AS1400'
            }
        
        if platform == 'Illumina':
            kit_details['DNEASY_KIT'] = {
                'title': 'DNeasy Kit',
                'vendor': 'Qiagen',
                'version': '',
                'catalog_number': '69504'
            }
        else:  # PacBio
            kit_details['MAGATTRACT_HMW_DNA_KIT'] = {
                'title': 'MagAttract HMW DNA Kit',
                'vendor': 'Qiagen',
                'version': '',
                'catalog_number': 67563
            }
            # Add PROMEGA kit for PacBio if BLOOD tissue exists
            if 'BLOOD' in tissue_types:
                kit_details['PROMEGA_RSC_BLOOD_DNA_KIT'] = {
                    'title': 'Promega RSC Blood DNA Kit',
                    'vendor': 'Promega',
                    'version': '',
                    'catalog_number': 'AS1400'
                }
    
    # Always include library prep kit
    if platform == 'Illumina':
        kit_details['KAPA-HYPERPREP-KIT'] = {
            'title': 'KAPA HyperPrep Kit',
            'vendor': 'Roche Diagnostics Systems',
            'version': '',
            'catalog_number': '7962371001'
        }
    if platform == 'RNAseq':
        kit_details['POLARIS-DEPLETION-KIT'] = {
            'title': 'Polaris Depletion Kit',
            'vendor': 'Watchmaker Genomics',
            'version': '',
            'catalog_number': ''
        }
    
    kit_records = []
    for kit_name, details in kit_details.items():
        kit_records.append({
            'submitted_id': format_preparation_kit_id(prefix, kit_name),
            'title': details['title'],
            'vendor': details['vendor'],
            'version': details['version'],
            'catalog_number': details['catalog_number']
        })
    
    return pd.DataFrame(kit_records)


def generate_treatment_sheet():
    """Generate treatment metadata sheet (empty for this dataset)."""
    return pd.DataFrame(columns=[
        'submitted_id', 'agent', 'description', 'concentration', 
        'concentration_units', 'duration', 'temperature'  # 250730- Reordered
    ])


def generate_library_sheet(df, analyte_sheet, prefix="WASHU", platform="Illumina"):
    """Generate library metadata sheet from sequencing data."""
    library_df = get_unique_libraries(df)
    
    analyte_map = dict(zip(analyte_sheet['external_id'], analyte_sheet['submitted_id']))
    
    assay_map = {
        'RNAseq': 'bulk_rna_seq',
        'Illumina': 'bulk_wgs_pcr_free',
        'PacBio': 'bulk_wgs_pcr_free'
    }
    assay = assay_map.get(platform, 'bulk_wgs_pcr_free')
    
    library_records = []
    for _, row in library_df.iterrows():
        sample_name = row['Sample Name']
        donor_id, tissue_code, sample_suffix = parse_sample_name(sample_name)
 
        analyte_id = ''
        for analyte_external_id, analyte_submitted_id in analyte_map.items():
            analyte_donor, analyte_tissue, _ = parse_sample_name(analyte_external_id)
            if analyte_donor == donor_id and analyte_tissue == tissue_code:
                analyte_id = analyte_submitted_id
                break

        # Add all required columns
        library_records.append({
            'submitted_id': format_library_id(prefix, sample_name, row['ESP ID']),
            'comments': '',
            'external_id': '',
            'description': '',
            'a260_a280_ratio': '',
            'adapter_name': '',
            'adapter_sequence': '',
            'amplification_cycles': '',
            'analyte_weight': '',
            'antibody': '',
            'barcode_sequences': '',
            'concatenated_reads': '',
            'fragment_mean_length': '',
            'guide_sequence': '',
            'insert_coefficient_of_variation': '',
            'insert_maximum_length': '',
            'insert_mean_length': '',
            'insert_minimum_length': '',
            'preparation_date': '',
            'dna_target': '',
            'target_fragment_size': '',
            'target_insert_maximum_length': '',
            'target_insert_minimum_length': '',
            'target_monomer_size': '',
            'analytes': analyte_id,
            'assay': assay,
            'library_preparation': get_library_preparation_id(prefix, platform)
        })
    
    return pd.DataFrame(library_records)


def generate_library_preparation_sheet(prefix="WASHU", platform="Illumina"):
    """Generate library preparation metadata sheet."""
    
    if platform == 'RNAseq':
        prep_id = f"{prefix}_LIBRARY-PREPARATION_WATCHMAKER_BULK-RNA-SEQ"
        return pd.DataFrame([{
            'submitted_id': prep_id,
            'description': '',
            'adapter_inclusion_method': '',
            'amplification_method': '',
            'fragmentation_method': '',
            'insert_selection_method': '',
            'enzymes': '',
            'rna_seq_protocol': 'Watchmaker',
            'size_selection_method': '',
            'strand': 'First Stranded',
            'trim_adapter_sequence': '',
            'preparation_kits': f"{prefix}_PREPARATION-KIT_POLARIS-DEPLETION-KIT",
            'treatments': ''
        }])
    elif platform == 'PacBio':
        prep_id = f"{prefix}_LIBRARY-PREPARATION_PACBIO_BULK-WGS"
        return pd.DataFrame([{
            'submitted_id': prep_id,
            'description': '',
            'adapter_inclusion_method': 'Ligation',
            'amplification_method': '',
            'fragmentation_method': '',
            'insert_selection_method': '',
            'enzymes': '',
            'rna_seq_protocol': '',
            'size_selection_method': '',
            'strand': '',
            'trim_adapter_sequence': '',
            'preparation_kits': prep_id,
            'treatments': ''
        }])
    else:  # Illumina
        prep_id = f"{prefix}_LIBRARY-PREPARATION_ILLUMINA_BULK-WGS"
        return pd.DataFrame([{
            'submitted_id': prep_id,
            'description': '',
            'adapter_inclusion_method': 'Ligation',
            'amplification_method': '',
            'fragmentation_method': '',
            'insert_selection_method': '',
            'enzymes': '',
            'rna_seq_protocol': '',
            'size_selection_method': '',
            'strand': '',
            'trim_adapter_sequence': '',
            'preparation_kits': f"{prefix}_PREPARATION-KIT_KAPA-HYPERPREP-KIT",
            'treatments': ''
        }])


def generate_sequencing_sheet(prefix="WASHU", platform="Illumina"):
    """Generate sequencing metadata sheet."""
    
    if platform == 'RNAseq':
        seq_id = f"{prefix}_SEQUENCING_NOVASEQXPLUS_150BP_100M"
        return pd.DataFrame([{
            'submitted_id': seq_id,
            'read_type': 'Paired-end',
            'target_read_length': 150,
            'additional_notes': '',  # NEW: Added additional_notes
            'flow_cell': '25B',
            'movie_length': '',
            'on_target_rate': '',
            'target_coverage': '',
            'target_read_count': 100000000,
            'target_monomer_length': '',
            'sequencer': 'illumina_novaseq_x_plus',
            'preparation_kits': ''
        }])
    elif platform == 'PacBio':
        seq_id = f"{prefix}_SEQUENCING_PACBIO_60X"
        return pd.DataFrame([{
            'submitted_id': seq_id,
            'read_type': 'Single-end',
            'target_read_length': 20000,
            'additional_notes': '',  # NEW: Added additional_notes
            'flow_cell': '',
            'movie_length': 30,
            'on_target_rate': '',
            'target_coverage': 60,
            'target_read_count': '',
            'target_monomer_length': '',
            'sequencer': 'pacbio_revio_hifi',
            'preparation_kits': f"{prefix}_PREPARATION-KIT_REVIO-SEQUENCING-PLATE"
        }])
    else:  # Illumina
        seq_id = f"{prefix}_SEQUENCING_NOVASEQX-150X"
        return pd.DataFrame([{
            'submitted_id': seq_id,
            'read_type': 'Paired-end',
            'target_read_length': 150,
            'additional_notes': '',  # NEW: Added additional_notes
            'flow_cell': '25B',
            'movie_length': '',
            'on_target_rate': '',
            'target_coverage': 150,
            'target_read_count': '',
            'target_monomer_length': '',
            'sequencer': 'illumina_novaseq_x_plus',
            'preparation_kits': ''
        }])

# ============================================================================
# FILE-BASED SHEETS
# ============================================================================

def generate_fileset_sheet(df, library_sheet, prefix="WASHU", platform="Illumina"):
    """Generate fileset metadata sheet from sequencing data."""
    library_df = get_unique_libraries(df)
    
    seq_id = get_sequencing_id(prefix, platform)
    
    fileset_records = []
    for _, row in library_df.iterrows():
        library_id = format_library_id(prefix, row['Sample Name'], row['ESP ID'])
        
        fileset_records.append({
            'submitted_id': format_fileset_id(prefix, row['Sample Name'], row['ESP ID']),
            'description': '',
            'submitter_comments': '',  # 250730 - NEW: Added submitter_comments
            'libraries': library_id,
            'sequencing': seq_id,
            'samples': ''
        })
    
    return pd.DataFrame(fileset_records)


def generate_unaligned_reads_sheet(df, fileset_sheet, prefix="WASHU", platform="Illumina"):
    """Generate unaligned reads metadata sheet from sequencing data."""
    file_df = get_file_data(df)
    
    # Create mapping from ESP ID to fileset submitted_id
    esp_to_fileset = {}
    for _, row in fileset_sheet.iterrows():
        fileset_id = row['submitted_id']
        esp_id = fileset_id.split('_')[-1]  # Gets "LIB072724-DIL01" from "WASHU_FILE-SET_SMHT009-3A-02_LIB072724-DIL01"
        esp_to_fileset[esp_id] = fileset_id
    
    unaligned_records = []
    for _, row in file_df.iterrows():
        filename = row['filename']
        esp_id_from_filename = filename.split('_')[0]  
        
        # Find the corresponding fileset
        fileset_id = esp_to_fileset.get(esp_id_from_filename, '')
        
        # Capitalize filename extension
        # filename_capitalized = filename.replace('.fastq.gz', '.fastq.gz')
        submitted_id = format_unaligned_reads_id(prefix, filename)
        
        # Generate actual filename with prefix for RNAseq (for filename column only)
        if platform == 'RNAseq':
            sample_name = row['sample_name']  
            esp_id = row['esp_id']  
            actual_filename = f"{sample_name}.{esp_id}.{filename}"
        else:
            actual_filename = filename

        # For R2 files, find their paired R1
        paired_with = ''
        if row['read_pair'] == 'R2':
            # Find corresponding R1 file
            r1_filename = filename.replace('_R2_', '_R1_')
            # r1_filename_cap = r1_filename.replace('.fastq.gz', '.fastq.gz')
            paired_with = format_unaligned_reads_id(prefix, r1_filename)
        
        # Add all required columns
        unaligned_records.append({
            'submitted_id': submitted_id,
            'filename': actual_filename,
            'submitted_md5sum': '',
            'data_category': 'Sequencing Reads',
            'data_type': 'Unaligned Reads',
            'n50': '',
            'flow_cell_barcode': '',
            'flow_cell_lane': '',
            'description': '',
            'read_pair_number': row['read_pair'],
            'file_format': 'fastq_gz',
            'file_sets': fileset_id,  # Now correctly maps to the fileset that ends with the same ESP ID
            'derived_from': '',
            'external_quality_metrics': '',
            'paired_with': paired_with,
            'software': ''
        })
    
    return pd.DataFrame(unaligned_records)


def generate_aligned_reads_sheet(df):
    """Generate aligned reads metadata sheet (empty for this dataset)."""
    return pd.DataFrame(columns=[
        'submitted_id', 'filename', 'submitted_md5sum', 'data_category', 'data_type', 'n50',
        'flow_cell_barcode', 'flow_cell_lane', 'description', 'alignment_details',
        'file_format', 'file_sets', 'reference_genome', 'derived_from',
        'external_quality_metrics', 'software'
    ])


def generate_variant_calls_sheet():
    """Generate variant calls metadata sheet (empty for this dataset)."""
    return pd.DataFrame(columns=[
        'submitted_id', 'data_category', 'data_type', 'filename', 'submitted_md5sum',
        'description', 'comparator_description', 'external_databases', 'filtering_methods',
        'mode', 'file_format', 'file_sets', 'reference_genome', 'derived_from',
        'external_quality_metrics', 'software'
    ])


def generate_software_sheet():
    """Generate software metadata sheet (empty for this dataset).""" 
    return pd.DataFrame(columns=[
        'submitted_id', 'category', 'title', 'version', 'description', 
        'binary_url', 'commit', 'source_url', 'gpu_architecture',  # 250730- CHANGED: gpu → gpu_architecture
        'model', 'modification_tags'
    ])

# 250730 - added new sheet
def generate_supplementary_file_sheet():
    """Generate supplementary file metadata sheet (empty for this dataset)."""
    return pd.DataFrame(columns=[
        'submitted_id', 'data_category', 'data_type', 'filename', 'submitted_md5sum',
        'description', 'haplotype', 'target_assembly', 'source_assembly', 
        'file_format', 'file_sets', 'software', 'derived_from', 
        'donor_specific_assembly', 'external_quality_metrics', 'reference_genome'
    ])
