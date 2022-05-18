# XEM_model
XEM_model is a Fortran script used to generate inclusive scattering cross-sections from nuclei.  The inelastic cross-section is produced using the F1F2## model of world data produced by Eric Christy.  The quasi-elastic contribution is produced using a Y-scaling fit to 6GeV Hall C X>1 data.  

There are two main uses for this code:
```
1.) Create Standalone Cross-Section Tables
2.) Inject the XEM_model into rc_externals via the SIGMODEL_CALC routine
```

# Initial Setup
You'll need to add a few directories in the parent directory including:
```
input
output
output/figures
```
To compile on a JLab computer, simply type `make` in the XEM_model directory.  This will create an executable called XEM_model, which is called by `run_XEM_model` and `run_many_XEM_model`.  

# Usage 
To create a standalone cross-section table, one must first generate an input table.  This can be done using the jupyter notebook `MakeTbl.ipynb`.  Alternatively, you can create your own routine as long as the input file generated follows the format, `{:6.3f} {:6.3f} {:6.3f} {:d} {:d}`, which is: `E_beam (GeV) E_prime (GeV) Theta (degree) A Z`.

The output file has the following columns:
```
y - y-scaling variable
A - Atomic Mass
Z - Atomic Number
Theta - Angle from input file (degrees)
E_prime - Scattered electron energy (GeV)
x - Bjorken Scaling variable
inelastic cross-section - Calculated from smearing the F1F2## model for nuclei (nb/MeV/sr)
quasi-elastic cross-section - Calculated using y-scaling (nb/MeV/sr)
```

The inclusion of the XEM_model in rc_externals is straight-forward, but requires more details to set up the rc-externals makefile.  This discussion is left for in person discussion, but generally rc-externals makes use of SIGMODEL_CALC.  One has to be careful with the fortran REAL and REAL*8 types as they are not cast directly to the correct type in most cmpilers.  