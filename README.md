# Deuteros

## An updated Deuteros 2.0 is available and can be found at https://github.com/andymlau/Deuteros_2.0. 

*Readme last updated 29/01/2020*

![alt text](https://github.com/andymlau/Deuteros/blob/master/readme_imgs/GUI_screenshot.png)

Deuteros is software for the rapid analysis and visualization of data from differential hydrogen deuterium exchange-mass spectrometry. Hydrogen deuterium exchange-mass spectrometry (HDX-MS) has emerged as a powerful technique for interrogating the conformational dynamics of proteins and their complexes. Currently, analysis of HDX-MS data remains a laborious procedure, mainly due to the lack of streamlined software to process the large datasets. We present Deuteros which is a standalone soft-ware designed to be coupled with Waters DynamX HDX data analysis software, allowing the rapid analysis and visualization of data from differential HDX-MS.

Reference: https://academic.oup.com/bioinformatics/article/35/17/3171/5288775

#### Requirements for installation:
- MATLAB Runtime Library (2018 and onwards)
- Ghostscript (see below)
- Windows/MacOS system

##### Ghostscript

Deuteros is available for any machine capable of running MATLAB, and requires **Ghostscript** for figure exporting.

Ghostscript for Windows: https://www.ghostscript.com/download/gsdnld.html
Ghostscript for Mac: with homebrew installed, use `brew install ghostscript`

##### MATLAB RunTime Library

The MATLAB RunTime library can be downloaded from the MathWorks website free of charge: https://uk.mathworks.com/products/compiler/matlab-runtime.html

---

### FAQ

#### 1. I want to use Deuteros but I don't have a license for MATLAB

Running Deuteros **does not** require a MATLAB license. Deuteros needs the MATLAB Runtime library (2018b and onwards) to be installed. The library can be downloaded and installed free of charge from MathWorks: https://uk.mathworks.com/products/compiler/matlab-runtime.html

#### 2. Deuteros reads my files successfully but when I try to export files, nothing happens? 

By default files are exported to the same directory as your input files, but export may fail silently when Deuteros doesn't have write permissions to these directories, e.g. a network or external drive. Try placing the input files into the Desktop for example and try again.  

#### 3. Generating the 'difference' file for Deuteros
1. To export difference data you should open a butterfly plot in DynamX.
2. Right click on the butterfly plot and click properties.
3. Ensure you have “difference index” selected (rather than relative uptake) and click ok.
4. In the top left, the upper dropdown should be your bound or mutant state, with the unbound or wildtype in the second dropdown.
5. Right click on the plot again and click copy data.
6. Paste this into excel and save as “CSV (Comma delimited) (.csv)”, there are other .csv formats to choose from but it must be this one for it to work as intended.

#### 4. I tried to export an image but nothing happened?

Ensure that **Ghostscript** is installed on your system (see below) and try again. Also see Q2 above. 

---

### History
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

---

## Usage

The easiest method of running Deuteros is to open the Deuteros.m script using MATLAB and press 'Run'. This will initialise the Deuteros GUI directly. If running without MATLAB, install Deuteros as a standalone program using one of the installers for Windows or Mac and following the installer instructions. 

Test files have been provided on this repository: 

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

The below tutorial video covers how to use Deuteros:

[![Deuteros Tutorial](http://img.youtube.com/vi/4DHuDrj2MPI/0.jpg)](http://www.youtube.com/watch?v=4DHuDrj2MPI "Deuteros Tutorial")

---

### Contact

Please email queries to andy.lau (at) kcl.ac.uk 

If reporting a bug, please include details of how the error occured, screenshots, as well as files needed to reproduce the error.
