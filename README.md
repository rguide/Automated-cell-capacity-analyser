% Automated Cell Capacity Analyzer
%
% Author: Roman Kula, MD
% Institution: Department of Physiology, Faculty of Medicine, Masaryk University,
% Brno, Czech Republic
%
% This script automatically and uniformly determines cell membrane capacitance across
% a set of cardiomyocytes. The input file should be in ABF (Axon Binary File) format,
% obtained from pCLAMP software.
%
% Details on determining cell membrane capacitance are provided in the Supplementary
% Material of Kula R, Bébarová M, Matejovič P, Šimurda J, Pásek M. (2020). 
% "Current density as a routine parameter for describing ionic membrane current:
% Is it always the best option?" Prog Biophys Mol Biol. 157:24-32.
% doi:10.1016/j.pbiomolbio.2019.11.011. Briefly, the input ABF file must contain 
% the transient current recorded after applying a subthreshold voltage step.
% The transient charge is then obtained by integrating the area under the 
% time-current curve. Corrections for residual current flow and voltage 
% drop across series resistance are applied based on Platzer, D., Zorn-Pauly, K. (2016).
% "Letter to the editor: Accurate cell capacitance determination from a single voltage
% step: A reminder to avoid unnecessary pitfalls." Am. J. Physiol. Heart Circ. Physiol.
% 311, H1072–H1073. https://doi.org/10.1152/ajpheart.00503.2016.
%
% This script was created in Matlab R2017a and requires the abfload function
% by Harald Hentschke (2024). abfload can be found on the MATLAB Central File
% Exchange: https://www.mathworks.com/matlabcentral/fileexchange/6190-abfload.

% Instructions:
% 1. In the current folder, ensure the following:
%    - Place `raw_data.mat` (filled according to the instructions in the README file).
%    - Place all `.abf` files you want to analyze.
%    - Create a new empty folder named "Graphs" for saving output figures.
%    - Save this MATLAB script in the same folder for optimal file path handling.

% 2. To improve the accuracy of the analysis, use a higher sampling
%    frequency (e.g., 100 kHz). The output graphs are designed to visually
%    estimate the goodness of fit. Use only the traces that can be
%    adequately fitted by a one-term exponential equation.

% 3. Example usage:
%    - To demonstrate how to correctly run the analyzer, download the
%      `raw_data.mat` file from GitHub along with the corresponding `.abf`
%      files (recorded on CHO cells to determine membrane capacitance from
%      subthreshold voltage steps).
%    - Place these files in the same folder as this script and launch the analyzer.

% Detailed Explanation of Each Column
% To run the Analyser, raw_data.mat file must be prepared according to the following instruction:

% 1. File Name (File_name) - Column 1
% Description: The name of the .abf file you wish to analyze.
% Format: String (including the file extension .abf).
% Example: 'cell1.abf', 'experiment_data.abf'

% 2. Sweeps to Average (Average_Sweeps_No) - Column 2
% Description: Specifies which sweeps (traces) from the .abf file to include in the averaging process.
% Format:
% Specific Sweeps: A numerical array of sweep numbers.
% Example: [1, 2, 3] (to average sweeps 1, 2, and 3)
% All Sweeps: The string 'all' to include all available sweeps.
% Example: 'all'

% 3. Gain (gain) - Column 3
% Description: The gain used during data acquisition. This is used to adjust the signal amplitude appropriately.
% Format: Numeric value (usually a positive number).
% Example: 1, 10, 100

% 4. Sampling Frequency (Frequency) - Column 4
% Description: The sampling frequency at which the data was recorded, in kilohertz (kHz).
% Format: Numeric value (positive number).
% Example: 10 (for 10 kHz), 50, 100

% 5. Cursor1 (Cursor1) - Column 5
% Description: The starting time (in milliseconds) for the baseline adjustment interval.
% Format: Numeric value (time in milliseconds).
% Example: 5.0, 10.0

% 6. Cursor2 (Cursor2) - Column 6
% Description: The ending time (in milliseconds) for the baseline adjustment interval.
% Format: Numeric value (time in milliseconds).
% Example: 10.0, 15.0

% 7. Cursor3 (Cursor3) - Column 7
% Description: The time (in milliseconds) marking the start of the voltage impulse.
% Format: Numeric value (time in milliseconds).
% Example: 0.5, 1.0
