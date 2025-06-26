## 1- Depth distribution
Describes the code and data used to set up the model changing
particles depth distribution over time, by simulating distributions that fit
with a dataset of many depth distribution profiles.
<br/>
<br/>
The *Data* folder stores the initial dataset of profiles of Calanus densities (Calanus_Vertical_Distribution_Dataset_m3.csv) and other files created form this one by several R scripts. 
<br/>
<br/>
The R code *1_modifications.R* shapes the original dataset and create a new one for each copepod class. 
<br/>
<br/>
The R code *2_depth_profiles.R* creates figures of observed depth distribution profiles of each copepod class. 
<br/>
<br/>
The R code *3_optimisation.R* is the heart of the depth distribution model. For each copepod class, for each depth distribution profile, it fits 2 Orsnstein - Ulhenbek equations two simulates the profile as well as possible. At the end it stores the optimal parameters, and figures representing the observed and simulated profiles.
<br/>
<br/>
The R code *4_regression.R* correlates the optimal parameters with environmental variables.
