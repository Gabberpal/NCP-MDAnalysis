<div align="justify">

# Calculation of RMSD traces for NCP trajectory

### Theoretical notes

The processing of the MD trajectory of the NCP particle includes the following steps:

1) overlay all MD frames by superimposing them onto the reference coordinates 3lz0 via the secondary-structure Cα atoms
   from the histone core. At this step, we assemble the nucleosome particle that may appear divided at the boundaries due to
   periodic boundary condition in the MD simulation.


2) calculate RMSD using different sets of atoms:
- (i) Cα atoms that have been used to overlay the frames 
- (ii) N1 and N9 atoms from nDNA nucleobases
- (iii) sets (i) and (ii) combined
- (iv) N1 and N9 atoms from the inner turn of nDNA, nucleotides from -38 to 38
- (v) N1 and N9 atoms from the outer turn of nDNA, nucleotides from -72 to -39 and from 39 to 72
    
### Run scripts

The scripts for calculation of RMSD traces were assembeled into a pipeline using make utility. To process
your own trajectory, you need to:
1) copy the template [analysis_template](analysis_template)
2) specify the parameters in [analysis_template/common.mk](analysis_template/common.mk)
3) type make.

### Run tests

We provide the templates to analyze short 10-ns trajectories recorded with Amber and Gromacs packages. Of note, the
results are for demonstration purpose only and cannot be used for interpretation of H4 tail dynamics.

```code-block:: bash
   # run the script to process Amber MD trajectory 
   cd example/ 
   make
```

See the [example_output](example_output)

</div>
