Content:
========
This folder contains the Matlab routines constituting all models and data 
needed to reproduce the analysis and figures shown in the manuscript:
'Sample-based modeling reveals bidirectional interplay between cell cycle progression and extrinsic apoptosis' 


Requirements / System:
======================
The analysis has been performed using MATLAB2017a.
In Fig_Supplement_HCT the toolbox ‚IntiQuan IQM Tools‚ (Schmidt2006) was used.


- Dead_Survivor_green.csv/Dead_Survivor_white.csv contain experimental measurements of NCI-H460 cells that received TRAIL in S/G2/M and G1, respectively
- Dead_Survivor_green_HCT.csv/Dead_Survivor_white_HCT.csv contain experimental measurements of NCI-H460 cells that received TRAIL in S/G2/M and G1, respectively
- Untreated_CellCycle.csv contains phase lengths of untreated NCI-H460 cells
- Untreated_CellCycle_HCT.csv contains phase lengths of untreated HCT-116 cells - Function_Import files extract the experimental data



- executable files are named as Figures (Fig1, Fig2A,…)
- for parameter estimation: only exemplary, save_data files must be executed first


All executable files were adapted to cell line NCI-H460. For HCT-116, change the following:
- In ‚Function_import_Data…‘:  max. Experiment duration to 24.5 h
- In ‚Import_Cycle_Data‘: Untreated_CellCycle_HCT.csv
- In ‚Mean_Cycle_Start‘: 'Dead_Survivor_white_HCT.csv' and 'Dead_Survivor_green_HCT.csv'
