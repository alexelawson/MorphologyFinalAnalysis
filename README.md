Welcome to **Microglia Morphology Analysis**!

Included you can find Alex Lawson Data Analysis process (including raw files). 

**In [Morphoglia](./Morphoglia/) you can find all the data analysis from the morphoglia pipeline.**

Within the [Morphoglia/Data](./Morphoglia/Data) folder you can find:
- `Morphology.csv` which is the output from the Morphoglia Pipeline for the PVN.
- `Morphology_HDBSCAN_30_0.1_150_5.csv` which is the output from the UMAP and HDBSCAN clustering on all data in the PVN.
- [UMAP_HDBSCAN](./Morphoglia/Data/UMAP_HDBSCAN) has the UMAP and HDBSCAN umaps/clustering from a couple of the trialed parameters.
- [Feature_Selection](./Morphoglia/Data/Feature_Selection) contains the output from the Random Feature Selection of Morphoglia
- [Plots](./Morphoglia/Data/Plots) contains the plots, regenerated in python to illustrate results (look at [Image Processing Scripts](./Image%20Processing%20Scripts/corr-plots-step-by-step.ipynb) for code to generate them. 

**In [Working Data Analysis](./Working%20Data%20Analysis/) you can find Alex Lawson's raw working data analysis from all the Ciernia Lab protocol.**

This includes both the unedited raw files as well as the final analysis/methods to produce stats used in the paper. 
- `final-stats-process.csv` is the final process
- `analysis01-06` are the process and includes the intitial processing of the data, creation of cluster percentages, as well as intital stats and graphing

**In [Data Files](./Data%20Files) you can find the data files used in the Ciernia Lab morphology analysis.** 

This includes both a data frame with the raw individual feature values, cluster, and first 2 PCAs as well as cluster percentage results/data. 
- `final-dataframe.csv` is the data frame/data that was used to create the cluster percentage data

**In [Image Processing Scripts](./Image%20Processing%20Scripts/) you can find the scripts used in fiji, python, and R to process images as well as create some of the output graphs for both the Ciernia Lab and Morphoglia analysis.**

The file `Edited-Binary-Conversion` is the script that was used to take the RGB images and threshold them. This was adapted Directly from the Ciernia Lab's image processing script with a few modifications. You can adjust the threshold or use Auto/Local thresholding instead if it better suits your image. Once this is run on your original images you can manually pre process or directly plug into the Ciernia Lab's Single Cell and Skeleton pipeline. We ran the 'Close-' function again on our single cell images after pre-processing to remove any atrifacts from processing, and close any insignificant holes as these have a huge impact on the skeletonization. 


