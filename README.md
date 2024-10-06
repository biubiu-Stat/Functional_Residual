# Reproducibility Materials
Code to reproduce simulation results, figures and real data analysis results from the paper "A unified framework for residual diagnostics in generalized linear models and beyond"

## Organization

### Illustration

- "Sim_Illustration.R" contains R code for Figure 1-3 in Section 2.3

### Simulation-code

#### Simulation for ordinal data

- "Sim_Ordinal_Data.R" contains R code for 

   - Figure 4-5 in Section 3.1;
   - Figure S1-S4 in Supplementary Materials B;
   - Figure S21 and S22 in Supplementary Materials F.

#### Simulation for count data

- "Sim_Count_Data.R" contains R code for 
   - Figure S5-S9 in Supplementary Materials B;
   - Figure S17-S19 in Supplementary Materials D.
   
#### Simulation Graphical Test

- "Sim_Graphical_test.R" contains R code for Figure S23 in Supplementary Materials G;
- "Sim_Power_Comparison.R" contains R code for producing the results displayed in Table S4 in Supplementary Materials G.

### Simulation-results

Contains the results from running the code in the simulation-code folder as described above. 


### Real-data

* "winequality-white.csv" contains the real data set in Section 5.1. Detailed information about this dataset can be accessed [here](https://archive.ics.uci.edu/dataset/186/wine+quality). 

* "hour.csv" contains the real data set in Section 5.2. Detailed information about this dataset can be accessed [here](https://archive.ics.uci.edu/dataset/275/bike+sharing+dataset). 

### Real-data-code
- "Whitewine.R" contains the code to run real data analysis in Section 5.1. It produces 
   - Figure S10-S12 in Supplementary Materials B;
   - Table S1 in Supplementary Materials C.


- "Bike.R" contains the code to run real data analysis in Section 5.2. It produces    
   - Figure S13-S16 in Supplementary Materials B;
   - Table S3 in Supplementary Materials C.


- "Bike_TimeSeries.R" contains the code to run the temporal effects analysis and produces Figure S20 in Supplementary Materials E. 


- "Robustness.R" contains the code to run the robustness analysis and produces Figure S24 in Supplementary Materials E. 

