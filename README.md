# GMXAnalysis
For the purpose of post progressing of MD carried by gromacs  
This repository is used to store some useful scripts for MD analysis.



# Content:
- DCCM.py: a python script to read xtc and generate the DCCM(Dynamical Cross-Correlation Matrix) automatically
- tpr2gro.py: generate the gro file from tpr file (tpr should be converted )  
- genxyz_gui.py: generate a set of gaussian input files for rigid scan  


# Usage:
## DCCM.py
> required packages:
> - pandas
> - matplotlib
> - numpy
> - sys
> - mdanalysis   
 
**command:**   
`python DCCM.py md.gro md.xtc`

## tpr2gro.py
you should use `gmx dump` command to convert `.tpr` file to readable `out` file.  
Then use the command below:  
**command:**    
`python tpr2gro.py md.out`

## genxyz_gui.py
This script is used to generate a set of gaussian input files from the molecule initial position to the end position for the purpose of rigid scanning.  no-gui version is also provided.

> required packages:
> - numpy
> - pandas
> - math

**command:**    
`python genxyz_gui.xyz`





