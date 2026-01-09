# generate-smaht-metadata-washu

This repo contains easy-to-use scripts that automatically generate comprehensive metadata spreadsheets for genomics data submission from raw sequencing center CSV files.


## Requirements
- Python 3.10.11 with
  - `pandas` 1.5.3
  - `openpyxl` 3.0.10
- Raw sequencing data information in CSV format

## Usage
```bash
python generate_metadata.py --work-dir [FOLDER_PATH] --in-csv [CSV_FILENAME] --platform [PLATFORM]
```

### Parameters
- `--work-dir`: Folder containing your CSV file (output saved here too)
- `--in-csv`: CSV filename 
- `--platform`: Choose `Illumina`, `RNAseq`, or `PacBio`

### Examples
```bash
# Illumina data
python generate_metadata.py --work-dir /path/to/SMHT001 --in-csv SMHT001_Illumina_data.csv --platform Illumina

# RNA sequencing
python generate_metadata.py --work-dir /path/to/SMHT001 --in-csv SMHT001_RNA_data.csv --platform RNAseq

# PacBio data
python generate_metadata.py --work-dir /path/to/SMHT001 --in-csv SMHT001_PacBio_data.csv --platform PacBio
```

## Output
Running the scripts will create the Excel workbook: `{work_dir}/{SAMPLE_NAME}_{PLATFORM}_{YYYYMMDD}.xlsx`

Example: `SMHT001_Illumina_20251030.xlsx`

That contains multiple metadata sheets required for DAC genomics data submission.
