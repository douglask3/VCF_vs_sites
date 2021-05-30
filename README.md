# Re-assess the MODIS Vegetation continuous Fields tree cover in tropical forest-savanna region.

Contains scripts for Adzhar et al. submission to Biogeosciences.

## Getting data
Some data is already in the repo. MODIS VCF data needs downloading. To get MODIS data fro TROBIT sites, run in R:
>> source(download_VCF_data.r)

To download MODIS VCF tiles, run in command line:
>> sh grab_VCF.sh

## Reproducing work

VCF vs TROBIT optimization is performed in the same script as the scatter plot, Fig. 1 and correction Figure, Fig. A4, In R, run:
>> source("scatter_machine.r")


Figures are saved in "figs/" directory. Parameter distribution samples are saved in "outputs".

To draw the map (Fig 2), run in R:
>> source("plot_gridded_experiments.r")

This also make Fig. A5


Barplots (FIg. 3( are made running:
>> source("reconstruct_curves.r")

(I'll come up with better names at some point)


Figure A2 is made using:
>> source("overlap.r")

while A3 using:

>> source("clumping.r")
 
