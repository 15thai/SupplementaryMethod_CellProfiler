# Supplementary Method

This is the supplementary method for article: "arXiv"

## Getting Started

Original pipeline could be found at [1], modified version that we used can be found at: github 
Python codes to extract the feature data automatically if not familiar with MySQL. Example of feature measurement is also provided to observation. 
Original dataset could be downloaded at https://data.broadinstitute.org/bbbc/BBBC025/,
Running the pipeline to extract the feature measurement:
-	Pipeline was tested on CellProfiler V2.1.1 (window 64-bits), downloaded at: https://cellprofiler.org/releases/.
-	Folder name of images are separated as “#_name” such as “38034_ERSyto”. The predicted/original images should follow the same structure. 


### Prerequisites
Python:
- matlibplot 3.0.2
- pandas 0.24.0
- numpy 1.16

### Steps

In the Input Modules:
-	Open CellProfiler => Open project => Choose: U2OS_pipeline_Project.cpproj
-	From Image module, select the folder(s) containing the illumination .mat files and then drag-and-drop the folder(s) into the File list panel, indicated by the words “Drop files and folders here”. The panel will update with the locations of the individual image files. 
-	From Image module, select the folder(s) containing the images (and then drag-and-drop the folder(s) into the File list panel, indicated by the words “Drop files and folders here”. The panel will update with the locations of the individual image files.
-	Check to see the corresponding illumination .mat files (based on # of folder) are loaded. Use group file filter as an option to process sub data if processing memory is limited.
-	From Metadata module, change “Metadata file location” to point to the location of “TAK004_TAK005_U2OS_96H.csv”.

Output setting:
-	Point the location to the output folder to store the results from the program (Pipeline-Output => View output settings). Otherwise, default is used.

Feature Measurement:
-	Click on “Analyze Images”, wait until it finishes.
-	Do similarly for both 5 original channel images and 2 input + 3 prediction images to create raw feature measurement storage in “U2OS_CPA_FM_5_OriginalChannel” and “U2OS_CPA_FM_2_and_3_Predict”, respectively. In each case, update the output setting according to both output folders. 
-	In each folders, remove all folders/files (boundary segmented images + *Experiment.csv) except the folders with name started with “Image”
-	Run **“auto_extract.py”** point to the location of the output folder to produce summary excel files with correct directory.

```
python auto_extract.py .\U2OS_CPA_FM_5_OriginalChannel .\U2OS_CPA_FM_5_OriginalChannel_Final
python auto_extract.py .\ U2OS_CPA_FM_2_and_3_Predict  .\U2OS_CPA_FM_2_and_3_Predict _Final
```
-	Run **“correlation_score.py”** to product Pearson correction checkerboard.
```

```
-	Run **“correlation_score_distribute.py”** to product Pearson correction checkerboard distributed in Cell, Nuclei, Cytoplasm.
```

```

## Built With

* [CellProfiler](https://cellprofiler.org/releases/) - Software

## Authors

* **Thanh Nguyen** - *Initial work* - [GitHub](https://github.com/32nguyen)
* **Anh Thai** - *Initial work* - [GitHub]()

## Acknowledgments

* This is the modified version from original pipeline found at [1]

## Supplementary References
[1] S. Singh, X. Wu, V. Ljosa, M. Bray, F. Piccioni, D. E. Root, J. G. Doench, J. S. Boehm, and A. E. Carpenter. "Morphological profiles of RNAi-induced gene knockdown are highly reproducible but dominated by seed effects." PLoS One 10, no. 7 (2015)

