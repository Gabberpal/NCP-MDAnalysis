## 1. change the path to SCRIPT_DIR according to your directory tree
SCRIPT_DIR:=../../../scripts/

## 2. specify MD traj parameters
# set trajectory path; for example md_setup/md_protocol/TIP4P-D_disp-replica-1/ff99SB-disp/wt/02_production/6_run/
TRAJECTORY_PATH:=../../../../short_example_trajectory/trajectory_AmberNetCDF/
# set path to trajectory reference pdb (of note, we need ); for example md_protocol/TIP4P-D_disp-replica-1/ff99SB-disp/wt/01_equil_histone_tails/0_prepare/wt.pqr
TRAJECTORY_REFERENCE_PATH:=../../../../short_example_trajectory/reference.pdb
# set trajectory length
TRAJECTORY_LENGTH:=10 # ns
# set step of printing out coordinates in trajectory
DT_NS=0.01 # time step of output for frame coordinates


## 3. specify fit parameters
FIT_LIMIT_NS=10 # ns
# you may specify logariphmic resampling of the correlation function (LAG_INDEX="log") and set the corresponding number of points from 0 ns to FIT_LIMIT_NS (N_LAG_POINTS).
# this step may prevent overfitting of the correlation function at large timescales
# If  you want to fit without logariphmic spacing, you should delete the these parameters (LAG_SPACING and N_LAG_POINTS)
LAG_SPACING="log"
N_LAG_POINTS=100


## 4. specify experimental parameters needed for calculation
NMR_FREQ=1200e6 # Hz
TUMBLING_TIME=163.4 # ns