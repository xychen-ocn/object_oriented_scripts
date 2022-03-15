# object_oriented_scripts
This is a repository with matlab scripts written in an object oriented way.

### Purpose: 

streamline two tasks (finding SST gradients and produce cloud image analysis)
with different measurements from ATOMIC.

### Status:

Ongoing, need to complete ATOMIC_GOES16 and ATOMIC_dataVisual

 - [x] added analysis with multiple matrix-pairs similar to the ones examined by Meroni et al. (2021?)
 - [ ] do analysis for cross-SST wind component for the disregarded RHB segments
 - [ ] understand temporal cloud fraction better, do a correlation analysis with the SST, or actually fluxes.


### To-be updated and check for sensitivity:
 - some of the selected RHB transects are not straight enough; I need to revist this and check if the results are sensitive to this issue. 
   - related parameters and functions in the **ATOMIC_dataProcess.m** are:

     - ```trajdir_cart``` in ```compute_distance_travelled```: revisit with central difference to calculate it;
     - selection criteria in ```select_straight_traj_segments```: revisit with computing the change of trajdir with central difference as well, I can also add a check at longer distance (5 points) to see if the direction change exceed threshold.

