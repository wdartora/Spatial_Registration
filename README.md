McIntire Lab - Spatial Omics Image Registration Pipeline

Overview
This repository contains scripts and documentation for a semi-automatic pipeline to register and integrate spatial omics data from mouse brain cross-sections. The workflow integrates spatial transcriptomics (Visium) and spatial metabolomics (DESI), using the Allen Mouse Brain Atlas as a spatial reference.

Workflow Steps
0. Manual Image Registration
Manually register the representative H&E image using:

Filebuilder.bat (QuickNII utility)

QuickNII

VisuAlign

1. Visium Data Preprocessing & Registration
Scripts/notebooks:

spatial registration/Visium_registration.ipynb

spatial registration/visium_st_seurat_updated.R

2. DESI Data Preprocessing & Registration
Scripts/notebooks:

spatial registration/DESI_registration.ipynb

spatial registration/DESI_processing.r

spatial registration/convert_desi_ion_to_matrix.ipynb

3. Spot Label Mapping
Spots2Labels.ipynb

or compact version (not recommended for execution):
spatial registration/Spots2Labels_compact_workflow_do_not_run.ipynb

Spatial Registration Pipeline
1. Image Preprocessing
Select a representative H&E image.

Resize it (~16 MP) per QuickNII guidelines.

Use Filebuilder.bat to format the image for QuickNII.

Align the Allen Mouse Brain Atlas to the H&E image in QuickNII.

Export .xml and load in VisuAlign.

Perform fine non-linear alignment and export atlas maps.

2. Visium Data Preprocessing
Load Visium data using Seurat.

If a slide contains two tissue sections, split accordingly.

Adjust high-res images and coordinate system for each section.

3. Image Landmarking
Use R to manually select anatomical landmarks in H&E and atlas images.

Recommended: <24 landmarks, mostly along the brain border and hippocampus.

Save as a coordinate table.

4. First Spatial Registration
Load landmarks and images using Python (cv2).

Use cv2.findHomography to align reference image to each brain section.

Overlay atlas to assign annotations and export annotated pixels.

5. Ion Data Preprocessing (DESI)
Convert raw DESI ion data to tabular format.

Identify metabolite features that resemble cortical structure.

Reconstruct tissue image using x/y coordinates and feature intensities.

6. Secondary Modality Alignment (DESI to Visium)
Identify landmarks in reconstructed DESI images.

Align the DESI coordinate system to the Visium space.

7. Second Spatial Registration
Register DESI images to Visium sections using Python (homography).

Map DESI feature data onto each Visium spot.

8. Integration
Merge transcriptomic (Visium) and metabolomic (DESI) data.

Perform normalization, quality control, and downstream multi-omics analysis.

Repository Structure
plaintext
Copy
Edit
/
├── Proj01/
│   ├── Integration of spatial modalities in brain tissue.pptx
│   └── VizDezI/
├── spatial registration/
│   ├── convert_desi_ion_to_matrix.ipynb
│   ├── DESI_processing.r
│   ├── DESI_registration.ipynb
│   ├── visium_st_seurat_updated.R
│   ├── Visium_registration.ipynb
│   └── Spots2Labels_compact_workflow_do_not_run.ipynb
├── DESI_processing.r
├── Spots2Labels.ipynb
├── ST_pipeline-DESI-Regions.ipynb
├── ST_pipeline-Visium.ipynb
├── visium_st_seurat.R
└── README.md
Requirements
R:
Seurat, ggplot2, dplyr, sp, etc.

Python:
opencv-python (cv2), numpy, pandas, matplotlib

Tools:

QuickNII

VisuAlign

Allen Mouse Brain Atlas reference files

Usage
Manual alignment of H&E images using QuickNII and VisuAlign.

Run Visium preprocessing using R and Python scripts in the spatial registration/ folder.

Run DESI preprocessing, reconstruct DESI images, and landmark.

Perform spatial alignment and export annotations for downstream multi-modal integration.

License
MIT License

Contact
For questions or collaborations, please open an issue or contact the McIntire Lab.

