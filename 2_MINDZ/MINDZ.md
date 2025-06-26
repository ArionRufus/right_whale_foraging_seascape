<br/>
## 2- MINDZ
Describes the code and data of the model *MINDZ* with *CIOPS-E* forcings. 
The FORTRAN model *MINDZ* uses currents forcings to simulate particles trajectories. 
This version also integrates the depth distribution model. Finally, it also 
includes R codes aiming at changing the outputs of MINDZ into 2D density data.
<br/>
<br/>
To run the model with 1 current input, you have to execute the file *MINDZ_execution.out* in the *run* folder. It will find the data in the *Data* folder, and execute the *MINDZ* model with the set of parameters written in *run.list.CIOPS-E* file. 
<br/>
The Pyton file *multiple_dates.py* was created to do run the MINDZ model several times, from a date to another. It will automatically change the *run.list.CIOPS-E* file and execute the model.
<br/>
<br/>
The scripts of the model are written in the *model* folder. The main script is *test_model_parallel.f90*, which call other scripts.
