"""
utils.py

Shared utilities and constants for metadata generation.
Contains common functions, ID formatters, and configuration constants.
"""

import pandas as pd

# ============================================================================
# CONSTANTS AND CONFIGURATION
# ============================================================================

# ID Templates
ID_TEMPLATES = {
    'donor': '{prefix}_DONOR_{donor_id}',
    'tissue': '{prefix}_TISSUE_{tissue_id}',
    'tissue_sample': '{prefix}_TISSUE-SAMPLE_{sample_id}',
    'analyte': '{prefix}_ANALYTE_{sample_id}_GDNA',
    'analyte_preparation': '{prefix}_ANALYTE-PREPARATION_{kit_name}',
    'preparation_kit': '{prefix}_PREPARATION-KIT_{kit_name}',
    'library': '{prefix}_LIBRARY_{library_id}',
    'library_preparation': '{prefix}_LIBRARY-PREPARATION_{prep_type}',
    'sequencing': '{prefix}_SEQUENCING_{seq_type}',
    'fileset': '{prefix}_FILE-SET_{library_id}',
    'unaligned_reads': '{prefix}_UNALIGNED-READS_{filename}',
    'aligned_reads': '{prefix}_ALIGNED-READS_{sample_id}_{segment}'
}

# 250730-Direct mapping from tissue codes to sample_sources tissue names
# Based on actual WashU metadata patterns
TISSUE_PATTERNS = {
    '3A': 'BLOOD',
    '3C': 'ESOPHAGUS', 
    '3E': 'ASCEN_COLON',
    '3G': 'DESCEN_COLON',
    '3I': 'LIVER',
    '3K': 'L_ADRENAL',
    '3M': 'R_ADRENAL',  # Inferred pattern
    '3O': 'AORTA',
    '3Q': 'LUNG',
    '3S': 'HEART',
    '3U': 'L_TESTIS',
    '3W': 'R_TESTIS',  # Inferred pattern
    '3Y': 'L_OVARY',
    '3AA': 'R_OVARY',  # Inferred pattern
    '3AD': 'UV_SKIN',
    '3AF': 'SKIN',
    '3AH': 'MUSCLE',
    '3AK': 'FRONTAL-LOBE',
    '3AL': 'TEMPORAL-LOBE',
    '3AM': 'CEREBELLUM',
    '3AN': 'L-HIPPOCAMPUS',
    '3AO': 'R-HIPPOCAMPUS'
}

# UBERON mappings for tissues (from correctOutput)
UBERON_MAP = {
    'BLOOD': 'UBERON:0013756',
    'UV_SKIN': 'UBERON:0004264',
    'SKIN': 'UBERON:0036149', 
    'MUSCLE': 'UBERON:0011907',
    'ESOPHAGUS': 'UBERON:0035216',
    'ASCEN_COLON': 'UBERON:0001156',
    'DESCEN_COLON': 'UBERON:0001158',
    'LIVER': 'UBERON:0001114',
    'LUNG': 'UBERON:0008952',
    'HEART': 'UBERON:0036285'
}

# Donor information mapping - update age and sex once we hear back
DONOR_INFO = {
    'SMHT006': {'age': '', 'sex': '', 'eligibility': 'Yes', 'hardy_scale': ''},
    'SMHT008': {'age': '', 'sex': '', 'eligibility': 'Yes', 'hardy_scale': ''},
    'SMHT009': {'age': 87, 'sex': 'Female', 'eligibility': 'Yes', 'hardy_scale': ''},
    'SMHT012': {'age': '', 'sex': '', 'eligibility': 'Yes', 'hardy_scale': ''},
    'SMHT015': {'age': '', 'sex': '', 'eligibility': 'Yes', 'hardy_scale': ''},
    'SMHT028': {'age': '', 'sex': '', 'eligibility': 'Yes', 'hardy_scale': ''},
    'SMHT029': {'age': '', 'sex': '', 'eligibility': 'Yes', 'hardy_scale': ''},
    'SMHT031': {'age': '', 'sex': '', 'eligibility': 'Yes', 'hardy_scale': ''},
    'SMHT039': {'age': '', 'sex': '', 'eligibility': 'Yes', 'hardy_scale': ''}
}

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def parse_sample_name(sample_name):
    """Parse sample name into components: donor_id, tissue_code, sample_suffix."""
    parts = sample_name.split('-')
    if len(parts) >= 3:
        donor_id = parts[0]
        tissue_code = parts[1]
        sample_suffix = parts[2]
        return donor_id, tissue_code, sample_suffix
    return None, None, None

def get_tissue_from_code(tissue_code):
    """Map tissue code directly to sample_sources tissue name."""
    return TISSUE_PATTERNS.get(tissue_code, 'UNKNOWN')

# 250730 - add new func for sample_sources query
def generate_sample_sources_id(donor_id, tissue_code, donor_prefix="NDRI"):
    """Generate sample_sources ID using direct tissue code mapping."""
    tissue_name = get_tissue_from_code(tissue_code)
    return f"{donor_prefix}_TISSUE_{donor_id}-{tissue_code}-{tissue_name}"

def get_uberon_id(tissue_name):
    """Get UBERON ID for tissue."""
    return UBERON_MAP.get(tissue_name, '')

def determine_category(tissue_name):
    """Determine sample category based on tissue type."""
    return "Liquid" if tissue_name == "BLOOD" else "Core"

def extract_library_id(esp_id):
    """Extract library identifier from ESP ID."""
    return esp_id

def get_kit_from_method(method):
    """Map library method to preparation kit."""
    if 'KAPA' in method.upper():
        return 'KAPA-HYPERPREP-KIT'
    return 'UNKNOWN-KIT'

def get_extraction_method_from_tissue(tissue_name):
    """Determine extraction method based on tissue type."""
    if tissue_name == "BLOOD":
        return "Magnetic Beads"
    return "Enzymatic"

# 250730 - updates
def get_preparation_kit_for_tissue(tissue_name, platform="Illumina"):
    """Get appropriate preparation kit for tissue type and platform."""
    if platform == 'RNAseq':
        return "RNEASY_KIT_PLUS"
    elif tissue_name == "BLOOD" and platform == 'Illumina':
        return "PROMEGA_RSC_BLOOD_DNA_KIT"
    elif tissue_name == "BLOOD" and platform == 'PacBio':
        return "PROMEGA_RSC_BLOOD_DNA_KIT"
    elif platform == 'Illumina':
        return "DNEASY_KIT"
    else:  # PacBio
        return "MAGATTRACT_HMW_DNA_KIT"

def get_donor_info(donor_id, field):
    """Get specific donor information field."""
    return DONOR_INFO.get(donor_id, {}).get(field, '')

def add_donor_info(donor_id, age=None, sex=None, eligibility=None, hardy_scale=None):
    """Add or update donor information."""
    if donor_id not in DONOR_INFO:
        DONOR_INFO[donor_id] = {}
    
    if age is not None:
        DONOR_INFO[donor_id]['age'] = age
    if sex is not None:
        DONOR_INFO[donor_id]['sex'] = sex
    if eligibility is not None:
        DONOR_INFO[donor_id]['eligibility'] = eligibility
    if hardy_scale is not None:
        DONOR_INFO[donor_id]['hardy_scale'] = hardy_scale

# ============================================================================
# ID FORMATTERS
# ============================================================================

def format_donor_id(prefix, donor_id):
    """Format donor submitted ID."""
    return f"{prefix}_DONOR_{donor_id}"

def format_tissue_id(prefix, donor_id, tissue_code, tissue_name):
    """Format tissue submitted ID."""
    return f"{prefix}_TISSUE_{donor_id}-{tissue_code}-{tissue_name}"

def format_tissue_sample_id(prefix, sample_name):
    """Format tissue sample submitted ID."""
    return f"{prefix}_TISSUE-SAMPLE_{sample_name}"

def format_analyte_preparation_id(prefix, kit_name):
    """Format analyte preparation submitted ID."""
    return f"{prefix}_ANALYTE-PREPARATION_{kit_name}"

def format_preparation_kit_id(prefix, kit_name):
    """Format preparation kit submitted ID."""
    return f"{prefix}_PREPARATION-KIT_{kit_name}"

def format_library_id(prefix, sample_name, library_esp_id):
    """Format library submitted ID."""
    return f"{prefix}_LIBRARY_{sample_name}_{library_esp_id}"

def get_library_preparation_id(prefix, platform):
    """Get library preparation ID based on platform."""
    if platform == 'RNAseq':
        return f"{prefix}_LIBRARY-PREPARATION_WATCHMAKER_BULK-RNA-SEQ"
    elif platform == 'PacBio':
        return f"{prefix}_LIBRARY-PREPARATION_PACBIO_BULK-WGS"
    else:  # Illumina
        return f"{prefix}_LIBRARY-PREPARATION_ILLUMINA_BULK-WGS"

def format_sequencing_id(prefix, seq_type):
    """Format sequencing submitted ID."""
    return f"{prefix}_SEQUENCING_{seq_type}"

def format_fileset_id(prefix, sample_name, library_esp_id):
    """Format fileset submitted ID."""
    return f"{prefix}_FILE-SET_{sample_name}_{library_esp_id}"

def format_unaligned_reads_id(prefix, filename):
    """Format unaligned reads submitted ID."""
    base_name = filename.replace('.fastq.gz', '.FASTQ.GZ')  # Capitalize extension
    return f"{prefix}_UNALIGNED-READS_{base_name}"

def format_aligned_reads_id(prefix, sample_id, segment):
    """Format aligned reads submitted ID."""
    return f"{prefix}_ALIGNED-READS_{sample_id}_{segment}"

def format_analyte_id(prefix, sample_name, suffix="_GDNA"):
    """Format analyte submitted ID with configurable suffix."""
    return f"{prefix}_ANALYTE_{sample_name}{suffix}"

def get_sequencing_id(prefix, platform):
    """Get sequencing ID based on platform."""
    if platform == 'RNAseq':
        return f"{prefix}_SEQUENCING_NOVASEQXPLUS_150BP_100M"
    elif platform == 'PacBio':
        return f"{prefix}_SEQUENCING_PACBIO_60X"
    else:  # Illumina
        return f"{prefix}_SEQUENCING_NOVASEQX-150X"


# ============================================================================
# DATA PROCESSING HELPERS
# ============================================================================

def load_sequencing_data(csv_file):
    """Load and process sequencing data from CSV file."""
    df = pd.read_csv(csv_file)
    df = df.dropna(subset=['Sample Name'])
    df['Sample Name'] = df['Sample Name'].astype(str)
    return df

def get_unique_samples(df):
    """Extract unique samples from sequencing data."""
    samples = df['Sample Name'].unique()
    sample_data = []
    
    for sample in samples:
        donor_id, tissue_code, sample_suffix = parse_sample_name(sample)
        if donor_id and tissue_code:
            tissue_name = get_tissue_from_code(tissue_code)
            sample_data.append({
                'sample_name': sample,
                'donor_id': donor_id,
                'tissue_code': tissue_code,
                'tissue_name': tissue_name,
                'sample_suffix': sample_suffix
            })
    
    return pd.DataFrame(sample_data)

def get_unique_libraries(df):
    """Extract unique libraries from sequencing data."""
    libraries = df.groupby(['ESP ID', 'Sample Name', 'Library Kit or Method']).first().reset_index()
    return libraries

def get_file_data(df):
    """Extract file information from sequencing data."""
    files = []
    for _, row in df.iterrows():
        if pd.notna(row['Read1 File']):
            files.append({
                'filename': row['Read1 File'],
                'esp_id': row['ESP ID'],
                'sample_name': row['Sample Name'],
                'read_pair': 'R1',
                'lane': row['Lane'],
                'flowcell': row['Flowcell ID']
            })
        if pd.notna(row['Read2 File']):
            files.append({
                'filename': row['Read2 File'],
                'esp_id': row['ESP ID'],
                'sample_name': row['Sample Name'],
                'read_pair': 'R2',
                'lane': row['Lane'],
                'flowcell': row['Flowcell ID']
            })
    return pd.DataFrame(files)