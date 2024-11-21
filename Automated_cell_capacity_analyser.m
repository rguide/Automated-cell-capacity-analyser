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
% Material of Kula R, Bébarová M, Matejoviè P, Šimurda J, Pásek M. (2020). 
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
% ---------------------------------------------------------------------------------------

% Initialize
close all;

% Define paths
data_path = 'C:\Users\Roman Kula\Documents\MATLAB';
save_path = fullfile(data_path, 'Graphs');

% Load data
load(fullfile(data_path, 'raw_data.mat'));

% Preallocate arrays based on the number of rows in raw_data
num_files = size(raw_data, 1);
Charge_Capacity = zeros(1, num_files); 
Charge_Value = zeros(1, num_files); 
Impulse_Start = zeros(1, num_files); 
Cell = cell(1, num_files); 
Value_I1 = zeros(1, num_files); 
Value_I01 = zeros(1, num_files); 
Value_I02 = zeros(1, num_files); 
Value_T1 = zeros(1, num_files); 
Value_Re = zeros(1, num_files); 
Value_Rm = zeros(1, num_files); 

for ii = 1:num_files
    % Extract data from raw_data
    File_name = raw_data{ii, 1}; 
    Average_Sweeps_No = raw_data{ii, 2}; 
    gain = raw_data{ii, 3};
    Frequency = raw_data{ii, 4}; 
    Cursor1 = raw_data{ii, 5}; 
    Cursor2 = raw_data{ii, 6}; 
    Cursor3 = raw_data{ii, 7}; 

    % Load data file
    file = abfload(fullfile(data_path, File_name));
    
    % Determine which sweeps to average
    if isequal(Average_Sweeps_No, 'all')
        sweeps_to_average = 1:size(file, 3);  % Use all available sweeps
    else
        sweeps_to_average = Average_Sweeps_No;
    end
    
    % Select and average sweeps
    selected_sweeps = file(:, :, sweeps_to_average);
    average_data = mean(selected_sweeps, 3);
    
    % Gain adjustment
    average_data(:, 1) = average_data(:, 1) / gain;
    
    % Voltage unit adjustment
    if abs(average_data(1, 2)) < 1
        average_data(:, 2) = average_data(:, 2) * 100;
    end

    % Calculate Cursor positions based on Frequency
    Cursor2 = round((Cursor2 * Frequency) + 1);
    Cursor1 = round((Cursor1 * Frequency) + 1);
    Cursor3 = round((Cursor3 * Frequency) + 1);
    U1 = average_data(Cursor3, 2);
    U2 = average_data(Cursor1, 2);
    
    % Adjust baseline
    baseline_table = table(average_data(Cursor3:Cursor2, 1));
    if U1 > U2
        p1 = find(abs(average_data(Cursor3:Cursor2, 2)) < (abs(U1) + 1));
    elseif U1 < U2
        p1 = find(abs(average_data(Cursor3:Cursor2, 2)) > (abs(U1) - 1));
    end
    Start_of_impulse = p1(end) + Cursor3;
    baseline_avg = mean(average_data(Cursor1:Cursor2, 1));
    data1 = table2array(baseline_table);
    data2 = data1 - baseline_avg;

    % Determine maxval and data_special
    if U1 > U2
        maxval = min(data1);
    elseif U1 < U2
        maxval = max(data1);
    end
    v = find(data1 == maxval); 
    data_special = data2(1:v(1, 1));
    
    % Vectorized perpendicular distance calculation
    xx = [1, size(data_special, 1)];
    yy = [data_special(1), data_special(end)];
    m = (yy(2) - yy(1)) / (xx(2) - xx(1));
    xc = (1:size(data_special, 1))';
    yc = data_special;
    yl = ((m .* xc) + (m^2 .* yc) - (m * xx(1)) + yy(1)) / (1 + m^2);
    xl = xc - m * (yl - yc);
    perpDist = (xl - xc).^2 + (yl - yc).^2;
    [val, idx] = max(perpDist);

    % Calculate integral and other values
    data3 = data2(idx:end);
    osa2 = ([0:length(data3)-1] / Frequency)';
    charge = cumtrapz(osa2, data3);
    if U1 > U2
        charge = min(charge);
    else
        charge = max(charge);
    end
    
    % Prepare data for fitting
    data_pro_fit = data1(v(1, 1):end);
    r = length(data_pro_fit);
    osa = ([0:r-1] ./ Frequency)';
    df1A = gradient(data_pro_fit, osa);
    if U1 > U2
        [~, idd] = max(df1A);
    elseif U1 < U2
        [~, idd] = min(df1A);
    end
    data_pro_fit = data_pro_fit(idd:end);
    r = length(data_pro_fit);
    osa = ([0:r-1] ./ Frequency)';

    % Curve fitting
    ftA = fittype('a*exp(-b*x)+c');
    optionA = fitoptions('Method', 'NonlinearLeastSquares', 'Robust', 'LAR', 'Algorithm', 'Levenberg-Marquardt', ...
                         'StartPoint', [0.9, 0.9, 0.9], 'MaxFunEvals', 6000000, 'MaxIter', 2000000);
    wA = fit(osa, data_pro_fit, ftA, optionA);
    coefsA = coeffvalues(wA);
   
    I01 = average_data(Start_of_impulse-1,1);  % Ensure Cursor3 - 1 is valid
    I02 = coefsA(3);
    T1 = 1 / coefsA(2);
    
    % Calculate corrected charge
    if U1 > U2
        charge = charge - (abs(I02 - I01) * T1);
    else
        charge = charge + (abs(I02 - I01) * T1);
    end

    % Calculate resistance and other parameters
    I1 = charge / T1; 
    Re = (U2 - U1) / (I1 + I02 - I01);
    km = (U2 - U1) / (Re * I1) - 1;
    Rm = Re / km;

    % Calculate delta and charge capacity
    if U1 > U2
        delta = (-(U1 - U2) / 1000) * (1 + (abs(I02 - I01) / I1));
    else
        delta = (-(U1 - U2) / 1000) * (1 - (abs(I02 - I01) / I1));
    end
    Charge_Capacity(ii) = charge / delta;

    % Plot and save figure
    aa = round((3.5*Frequency) + 1);
    data_pro_y1 = data1(1:aa);
    r = length(data_pro_y1);
    data_pro_x1 = ([0:r-1]./Frequency)';
    if U1 > U2
        maxvalue2 = min(data_pro_y1);
    elseif U1 < U2
        maxvalue2 = max(data_pro_y1);
    end
    [soux souy] = find(data_pro_y1 == maxvalue2);
    data_pro_x2 = ((soux-1):r)./Frequency;
    soux = (soux-1)/Frequency;
    osa_pro_fit = [-(idd-1):(length(data_pro_x2)-idd)]./Frequency;
    ffA = coefsA(1)*exp(-coefsA(2)*osa_pro_fit) + coefsA(3);
    
    figure;
    plot(data_pro_x1,data_pro_y1,data_pro_x2,ffA);
    title(File_name)
    xlabel('Time (ms)')
    ylabel('Current (nA)')
    h=num2str(ii);
    print(h,'-dpng')
    saveas(gcf, fullfile(save_path, sprintf('%d.png', ii)));
    close;

    % Store results in preallocated arrays
    Charge_Capacity(ii) = charge / delta;
    Charge_Value(ii) = charge;
    Impulse_Start(ii) = (Cursor3 - 1) / Frequency;
    Cell{ii} = File_name;
    Value_I1(ii) = I1;
    Value_I01(ii) = I01;
    Value_I02(ii) = I02;
    Value_T1(ii) = T1;
    Value_Re(ii) = Re;
    Value_Rm(ii) = Rm;
end

% Final results table
Cell = Cell(~cellfun('isempty', Cell));
Results = table(Cell', Charge_Capacity', Charge_Value', ...
                Value_Re', Value_Rm', Value_I1', Value_I01', Value_I02', ...
                Impulse_Start', Value_T1', ...
                'VariableNames', {'Cell', 'Charge_Capacity', 'Charge_Value', ...
                'Re', 'Rm', 'I1', 'I01', 'I02', 'Impulse_Start', 'T1'});
            
% Save Results to Excel File
output_file = fullfile(data_path, 'Results.xlsx');
writetable(Results, output_file);

% Confirm the file save
disp(['Results have been successfully saved to: ', output_file]);
