# Deuteros

#### History
###### Version 1.08 (10/12/2018)
- Added Mac and Windows installers for v1.08
- Fixed bug where last accessed folder would not be remembered
- Critical values 1 and 2 changed to dropdown menus - users directly select the percentage confidence and the corresponding t-test distribution table values (two-tailed, df=2) are selected automatically
- Plotting style of Flattened Data Map changed to use Matlab's polygon function rather than plotting lines with linewidths
- Plotting style of Woods plot now better utilise the white space in the Woods Plots section of the GUI
- Colour scheme of Woods plots updated
- Woods plot legends moved to inside to maximise the size of subplots
- The percentage confidence is now shown at the top left of each subplot
- Removed x-axis and y-axis from each subplot, there is one for the whole panel now
- Can now format more than one chain by separating chain IDs with commas in the PDB chain box (e.g. "A,B,C,D")
- Reproportioned the GUI to maximise the plotting space

###### Version 1.07
- Added number of peptides to Flattened Data Map
- Upgraded legends on Woods plot
- Fixed bug with importing data with modifications in difference file
- Fixed bug with exporting Woods plot without figure legends
- Fixed bug with y-axis limits for Woods plots

## Summary

![alt text](https://github.com/andymlau/Deuteros/blob/master/readme_imgs/GUI_screenshot.png)

Software for rapid analysis and visualization of data from differential hydrogen deuterium exchange-mass spectrometry

See also: https://www.biorxiv.org/content/early/2018/09/14/417998

Andy M. C. Lau, Zainab Ahdash, Chloe Martens, Argyris Politis

King's College London, Department of Chemistry, Britannia House, 7 Trinity Street, London, SE1 1DB, United Kingdom

Hydrogen deuterium exchange-mass spectrometry (HDX-MS) has emerged as a powerful technique for interrogating the conformational dynamics of proteins and their complexes. Currently, analysis of HDX-MS data remains a laborious procedure, mainly due to the lack of streamlined software to process the large datasets. We present Deuteros which is a standalone soft-ware designed to be coupled with Waters DynamX HDX data analysis software, allowing the rapid analysis and visualization of data from differential HDX-MS.

Deuteros is available for any machine capable of running MATLAB, and requires **Ghostscript** for figure exporting.

Ghostscript for Windows: https://www.ghostscript.com/download/gsdnld.html

Ghostscript for Mac: with homebrew installed, use `brew install ghostscript`

## Usage

The easiest method of running Deuteros is to open the Deuteros.m script using MATLAB and press 'Run'. This will initialise the Deuteros GUI directly. 

### Inputs
```
XylE WT apo vs E397Q difference.csv                     DynamX 'difference' file
XylE WT apo vs E397Q state.csv                          DynamX 'state' file
```

### Outputs
```
XylE WT apo vs E397Q difference_peptide_list.csv        List of deprotected and protected peptides  
XylE WT apo vs E397Q difference coverage.pdf            Vector: global coverage
XylE WT apo vs E397Q difference redundancy.pdf          Vector: peptide redundancy
XylE WT apo vs E397Q difference heatmap 5 min.pdf       Vector: heatmap of global deuterium uptake for t=5min
XylE WT apo vs E397Q difference_WoodsPlot.pdf           Vector: grid of Woods plots 
XylE WT apo vs E397Q difference_coverage.pml            PyMOL: coverage data for PDB
XylE WT apo vs E397Q difference_redundancy.pml          PyMOL: redundancy data for PDB 
XylE WT apo vs E397Q difference_uptake...scale.pml      PyMOL: Woods plot data (per timepoint) for PDB
```

### Tutorial Video

[![Deuteros Tutorial](http://img.youtube.com/vi/4DHuDrj2MPI/0.jpg)](http://www.youtube.com/watch?v=4DHuDrj2MPI "Deuteros Tutorial")
