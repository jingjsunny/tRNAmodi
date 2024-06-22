# tRNAmodi

# Project Title

Large-scale LCMS MRM data process workflow

## Description
This repository contains R scripts for processing and analyzing large-scale LCMS/MS MRM data. These scripts automate the data cleaning, evaluation, normalization, and fold change calculations necessary for high-throughput mass spectrometry data analysis.

## Getting Started

### Input Files
1. The input for these scripts should be a CSV datasheet originating from dynamic MRM or MRM data, pre-processed in Masshunter Quantitative Analysis software. Example datasets are provided here for reference. Modify the script paths to point to your data directory.
2. A CSV file titled "modification_inclusion_list.csv" is utilized for the calculation of fold changes. This file helps to filter out modifications that are unique to specific PA14 mutants, ensuring they are not included in the fold change analysis.
3. An Excel file (in ".xlsx" format) containing comprehensive details about mutant locations and PA14 gene annotations is required. It is integrated with the normalized peak area and fold change data output by the scripts.

### Data Preprocessing
Data preprocessing in Masshunter Quantitative Analysis software should be performed according to the specific method, such as peak selection and retention time alignment. Instructions for this preprocessing step are method-dependent and are not covered in this repository.

### Scripts and Functions
The scripts read data files, clean and process the LCMS data, and calculate normalized peak areas, fold changes, and detect outliers. Below is a brief overview of the main steps in the workflow: 
1. Read and prepare the batch table from CSV files. 
2. Clean up data and prepare column names. 
3. Extract MS and UV signals, and calculate normalized peak areas. 
4. Identify samples and peaks wrongly selected with retention times. 
5. Calculate fold change, filtering outliers based on defined thresholds.
6. Annotation of sample location and gene name.

### Output
The script will generate the csv files as following output:

1. Normalized peak areas for all samples ('Peakarea_all')
2. Fold change data for all samples ('FC_all')
3. Data for fold change above 2 or below 0.5 ('FC_F2_data')
4. List of peaks with varied retention time ('RT_ab_list')
5. List of samples with varied UV signals ('UV_ab_all')


