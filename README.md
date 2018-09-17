# Deuteros
Software for rapid analysis and visualization of data from differential hydrogen deuterium exchange-mass spectrometry

See also: https://www.biorxiv.org/content/early/2018/09/14/417998

Andy M. C. Lau, Zainab Ahdash, Chloe Martens, Argyris Politis

King's College London, Department of Chemistry, Britannia House, 7 Trinity Street, London, SE1 1DB, United Kingdom

![alt text](https://github.com/andymlau/Deuteros/blob/master/readme_imgs/GUI_screenshot.png?raw=true)

Hydrogen deuterium exchange-mass spectrometry (HDX-MS) has emerged as a powerful technique for interrogating the conformational dynamics of proteins and their complexes. Currently, analysis of HDX-MS data remains a laborious procedure, mainly due to the lack of streamlined software to process the large datasets. We present Deuteros which is a standalone soft-ware designed to be coupled with Waters DynamX HDX data analysis software, allowing the rapid analysis and visualization of data from differential HDX-MS.

Deuteros is available for any machine capable of running MATLAB, and requires **Ghostscript** for figure exporting.

Ghostscript for Windows: https://www.ghostscript.com/download/gsdnld.html

Ghostscript for Mac: with homebrew installed, use `brew install ghostscript`

## Usage

The easiest method of running Deuteros is to open the Deuteros.m script using MATLAB and press 'Run'. This will initialise the Deuteros GUI directly. 

### Inputs
```
XylE WT apo vs E397Q difference.csv                         DynamX 'difference' file
XylE WT apo vs E397Q state.csv                              DynamX 'state' file
```

### Outputs
```
XylE WT apo vs E397Q difference_peptide_list.csv            List of deprotected and protected peptides  
XylE WT apo vs E397Q difference coverage.pdf                Vector: global coverage
XylE WT apo vs E397Q difference redundancy.pdf              Vector: peptide redundancy
XylE WT apo vs E397Q difference heatmap 5 min.pdf           Vector: heatmap of global deuterium uptake for *t* = 5min
XylE WT apo vs E397Q difference_WoodsPlot.pdf               Vector: grid of Woods plots 
XylE WT apo vs E397Q difference_coverage.pml                PyMOL: coverage data for PDB
XylE WT apo vs E397Q difference_redundancy.pml              PyMOL: redundancy data for PDB 
XylE WT apo vs E397Q difference_uptake...scale.pml          PyMOL: Woods plot data (per timepoint) for PDB
```




Coming shortly!

