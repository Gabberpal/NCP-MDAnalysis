## 1. change the path to SCRIPT_DIR according to your directory tree
SCRIPT_DIR:={YOUR PATH}/NCP-MDAnalysis/RMSD/scripts/

## 2. specify MD traj parameters
# set trajectory path; for example md_setup/md_protocol/TIP4P-D_disp-replica-1/ff99SB-disp/wt/02_production/6_run/
TRAJECTORY_PATH:=""
# set path to trajectory reference pdb (of note, we need ); for example md_protocol/TIP4P-D_disp-replica-1/ff99SB-disp/wt/01_equil_histone_tails/0_prepare/wt.pqr
TRAJECTORY_REFERENCE_PATH:=""
# set trajectory length
TRAJECTORY_LENGTH:="" # ns
# set step of printing out coordinates in trajectory
DT_NS=0.01 # time step of output for frame coordinates
# set trajectory stride (DT_NS * TRAJECTORY_STRIDE is frequency of saving results)
TRAJECTORY_STRIDE:=100 # e.g. the result of RMSD calculation is saved every 1 ns (DT_NS=0.01 * TRAJECTORY_STRIDE:=100)

## 3. set path to Xray reference pdb relative to which RMSD is calculated
XRAY_REF="" # for example 3lz0