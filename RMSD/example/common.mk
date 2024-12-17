## 1. change the path to SCRIPT_DIR according to your directory tree
SCRIPT_DIR:=../../scripts/

## 2. specify MD traj parameters
# set trajectory path; for example md_setup/md_protocol/TIP4P-D_disp-replica-1/ff99SB-disp/wt/02_production/6_run/
TRAJECTORY_PATH:=../../../short_example_trajectory/trajectory_AmberNetCDF/
# set path to trajectory reference pdb (of note, we need ); for example md_protocol/TIP4P-D_disp-replica-1/ff99SB-disp/wt/01_equil_histone_tails/0_prepare/wt.pqr
TRAJECTORY_REFERENCE_PATH:=../../../short_example_trajectory/reference.pdb
# set trajectory length
TRAJECTORY_LENGTH:=10 # ns
# set type of trajectory files "nc" - AmberNetCDF, "xtc" - GromacsXtcFile
FILETYPE:=nc
# set pattern of trajectory files files: run00001.dat ---> "run%05d".dat
PATTERN:=run%05d
# set step of printing out coordinates in trajectory
DT_NS=0.01 # time step of output for frame coordinates
# set trajectory stride (DT_NS * TRAJECTORY_STRIDE is frequency of saving results)
TRAJECTORY_STRIDE:=100 # e.g. the result of RMSD calculation is saved every 1 ns (DT_NS=0.01 * TRAJECTORY_STRIDE:=100)

## 3. set path to Xray reference pdb relative to which RMSD is calculated
XRAY_REF=../../../short_example_trajectory/3lz0.pdb