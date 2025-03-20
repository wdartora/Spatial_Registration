# McIntire Lab - Spatial Omics Image Registration Pipeline

## Overview
This repository contains scripts and documentation for a semi-automatic pipeline to register and integrate spatial omics data from cross-sections of mouse brain. The workflow integrates spatial transcriptomics (*Visium*) and spatial metabolomics (*DESI*), using the *Allen Mouse Brain Atlas* as a reference.

## Workflow Steps

### 0. Manual Image Registration
Manually register the H&E image using:
- `Filebuilder.bat`
- `QuickNII`
- `VisuAlign`

### 1. Visium Data Preprocessing & Registration
Scripts:
- `01a_c_Visium_Data_Preprocessing_and_Post_Registration.R`
- `01b_Visium_Spots_Registration.py` (called in the middle of the R script)

### 2. DESI Data Preprocessing & Registration
Scripts:
- `02a_DESI_Preprocessing.py`
- `02b_d_DESI_Script.R`
- `02c_DESI.py` (called in the middle of the R script)

## Spatial Registration Pipeline

### 1. Image Preprocessing
- Parse and select a representative H&E image.
- Adjust the image size to fit *QuickNII* recommendations (~16 MP resolution).
- Prepare input files using `QuickNII`'s `Filebuilder.bat`.
- Align the 3D cross-section from *Allen Mouse Brain Atlas* to the H&E image.
- Export XML from *QuickNII* and use as input for *VisuAlign*.
- Perform non-linear transformations in *VisuAlign*.
- Export atlas maps with brain region annotations.

### 2. Visium Data Preprocessing
- Import spatial *Visium* data using *Seurat*.
- If necessary, split slides into two separate tissue sections.
- Adjust high-resolution images and spot coordinates to match coordinate system.

### 3. Image Landmarking
- Use R to interactively select landmarks on H&E and atlas images.
- Limit to ~24 landmarks, emphasizing brain borders and key structures (e.g., hippocampus).
- Convert landmarks into a tabular dataset.

### 4. First Spatial Registration
- Import landmarks and images in Python (`cv2`).
- Use *findHomography* (least squares) to align reference atlas to brain sections.
- Assign each pixel in the image section to the corresponding annotation region.
- Export pixel annotations as a tabular dataset.

### 5. Ion Data Preprocessing (DESI)
- Convert raw *DESI* data into a tabular format.
- Identify control features reflecting cortical structure.
- Convert DESI control features into images using x-y coordinates and intensity values.

### 6. Secondary Modality Alignment (DESI to Visium)
- Align *DESI* data to *Visium* coordinate space.
- Identify landmarks in *DESI* images of control features.

### 7. Second Spatial Registration
- Align *DESI* control feature images to *Visium* sections using Python.
- Link each *Visium* spot coordinate to the corresponding *DESI* data.
- Enable integration of transcriptomic and metabolomic data.

### 8. Integration
- Merge *DESI* and *Visium* data.
- Perform normalization and cleaning.
- Conduct downstream multiomic analysis.

## Repository Structure
```plaintext
/
├── scripts/
│   ├── 01a_c_Visium_Data_Preprocessing_and_Post_Registration.R
│   ├── 01b_Visium_Spots_Registration.py
│   ├── 02a_DESI_Preprocessing.py
│   ├── 02b_d_DESI_Script.R
│   ├── 02c_DESI.py
│   └── utilities/
├── data/
│   ├── raw/
│   ├── processed/
├── results/
├── README.md
└── LICENSE
```

## Requirements
- **R** (Seurat, ggplot2, dplyr, etc.)
- **Python** (cv2, numpy, pandas, etc.)
- **QuickNII** & **VisuAlign**
- *Allen Mouse Brain Atlas*

## Usage
1. Manually register H&E images (`QuickNII`, `VisuAlign`).
2. Run `01a_c_Visium_Data_Preprocessing_and_Post_Registration.R`.
3. Run `02a_DESI_Preprocessing.py` and subsequent steps.
4. Perform spatial alignment and integrate datasets.

## License
This project is licensed under the MIT License.

## Contact
For questions or collaborations, please open an issue or reach out via email.
