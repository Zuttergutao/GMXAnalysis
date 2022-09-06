# GMXAnalysis
For the purpose of post progressing of MD carried by gromacs  
This repository is used to store some useful scripts for MD analysis.

---

# Content:
- DCCM.py: a python script to read xtc and generate the DCCM(Dynamical Cross-Correlation Matrix) automatically
- tpr2gro.py: generate the gro file from tpr file (tpr should be converted )

---
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

## tpr2gro
you should use `gmx dump` command to output readable file to use.
