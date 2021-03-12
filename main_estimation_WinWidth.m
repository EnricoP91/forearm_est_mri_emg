close all
clear all
clc

%% SIMULATION PARAMETERS

ifshowplot = true;

% warning suppressing 
 id = 'MATLAB:gui:array:InvalidArrayShape';
 warning('off', id);

%% PATH SETUP
cur_fld = pwd;               
addpath(fullfile(cur_fld, 'lib'));
emg_path = fullfile(cur_fld, 'data', 'emg');        % filtered emg data
mri_path = fullfile(cur_fld, 'data', 'mri');        % information from segmentation    
graph_path = fullfile(cur_fld, 'data', 'graph');    % graph structures (incidence matrix, Lagrangian, etc)

%% Loading general model files
% *************************** INFORMATION *********************************
% 
% The file "est_voltage_active_DISTANCES" contain the matrices
%     B             : incidence matrix of the network
%     I_star_mat    : matrix of injected muscle currents (1 column = activation pattern) 
%     R             : diagonal matrix of resistances create using the distances measured with ImageJ 
%     R_inv         : inverse of resistance matrix R
%     LR            : Conductance matrix of the newtork 
%     G_mat         : shunt conductance matrix   

% The file "mat_electric_net_DISTANCES" contain the matrices
%     V                  : voltages measured on the nodes with single muscle activated one by one
%     electrode_index    : vector indicating the electrode nodes
%     muscle_index       : vector indicating the muscle nodes
% 
% *************************************************************************

% Area matrix loading (obtained from the graph creation and from the segmentation)
fprintf("Loading graph files . . .\n")
numbering = load(fullfile(graph_path, 'est_voltage_active_DISTANCES_AREA.mat'));  
graph_info = load(fullfile(graph_path, 'mat_electric_net_DISTANCES_AREA.mat'));
fprintf("Done.\n\n")

% Measurements of muscles from segmentation
fprintf("Loading muscles measurements. . .\n")
musclemeasurements = importCsvFile(fullfile(mri_path, 'measurements.csv'));
fprintf("Done.\n\n")

% EMG data
emg_file = uigetfile('*.mat', 'Select a EMG File', fullfile(emg_path, 'wext_semg.mat'));
fprintf("Loading emg files. . .\n")
emg_data = load(fullfile(emg_path, emg_file));
fprintf("Done.\n\n")

%% 0)                   GENERAL PARAMETERS
% Graph nodes 
n_node = size(numbering.muscle_index,1);
n_muscle = size(musclemeasurements,1)-1;
n_electrode = n_node-n_muscle;

% EMG information
semg_selected = emg_data.semg_Aprox_filt;
[n_channel, n_samples] = size(semg_selected);
    
%% 1)                       FORWARD PROBLEM
% =========================================================================
% Given input unit current matrix I_star_mat, Calculate the voltage profile
% on the electrodes for each of those nodes (muscles) active. This allow the
% calculation of the LFM which is SUBJECT-DEPENDANT. It describe the system
% impulsive response as it is described by the model
% =========================================================================     

% MRI electrode order list - ordering the electrodes according to anatomy
% (the electrodes indices around the arm are not in order, so they need to 
% be put in the right order)

selected_el = [1:n_electrode];       % order of electrodes along the arm which is different from the order in the graph

% Order of muscles. It needs to be the same with the graph one (the circumference of the arm is ignored)
muscles = categorical(musclemeasurements.Label);
muscles(muscles == "arm_circ") = [];

% Fixing eventual high conditioning values to avoid problem in in matrix inversion
fprintf('Checking the system matrix for high condition index value. \n')
[LFM,LFM_full] = fix_lfm(graph_info.LR, numbering.electrode_index);
fprintf('Checking complete. \n\n')

%% 2)                       INVERSE PROBLEM
% =========================================================================
% Given the ouput value of the voltage, we want to see whether we can 
% estimate the actives nodes, and if this estimation makes sense 
% =========================================================================
% @ input 
%           - the "lead field matrix" : formed with the voltages on single
%           activation of muscles on the columns (the order IS IMPORTANT)
%           - the read sEMG that will be fitted
% @ output 
%           - the weights at each muscles (the current at each muscle)
% -------------------------------------------------------------------------
% Definition of the windowing parameters
WINDOW_TIME_MS = 100;
OVERLAP_PERCENTAGE = 75;                % [%]
Fs = 2000;                              % EMG sampling frequency [Hz]
Ts = 1/Fs;                              % EMG sampling period  [sec]   
win_time = WINDOW_TIME_MS*1e-3;         % [sec]
% -------------------------------------------------------------------------
% EMG time vector
emg_data.time = [0 : Ts : (n_samples-1)*Ts];
% Creation of the window matrices        
width_win = round(win_time/Ts);
overlap = round(width_win * OVERLAP_PERCENTAGE/100);
step = width_win - overlap;
    
n_windows = round((n_samples - overlap)/step)-1;                % number of full windows 
lastwind_width = n_samples - (n_windows-1)*step - width_win;     % length of the last windows
window_start = [1 : step : (step*(n_windows))];                  % vector with the first sample of the window
window_end = window_start+(width_win-1);                          % vector with the last sample of the window
% adding last window
window_end = [window_end n_samples];
window_start = [window_start, (n_samples - lastwind_width - overlap)];
n_windows = n_windows + 1;
    
% Normalization of the cross-section of the muscles by the highest value
indtmp = 1;
for musID = 1 : height(musclemeasurements)  
    if musclemeasurements.Label(musID) ~= "arm_circ"
        indmus = find(muscles == musclemeasurements.Label(musID));
        fprintf("muscle : %s \n index : %d \n\n", musclemeasurements{musID,1},indmus)
        muscleareas(indmus,1) = musclemeasurements.Area(musID);
    end
end
muscleareas_norm = muscleareas./max(muscleareas);

% Figure set up
if ifshowplot == true
    figure_info = prepare_figure(semg_selected, [n_electrode, n_muscle], muscles, emg_data.time, win_time);
end
%% ESTIMATION LOOP
fprintf("Estimation over the windows \n")
VnodesMat =[];
IestMat = [];
VoutMat = [];
rms_semgRealMat = [];

for i = 1 : n_windows
        
         % window selection and rms calculation
         sEMG_window = semg_selected(:,(window_start(i):window_end(i)));
         rms_semgwindow = rms(sEMG_window,2);

         % rms with electrodes according to graph order (consider that there are 2 electrodes sheets)
         rms_semggraph = zeros(n_electrode,1);   
         for elecID = 1 : n_electrode
            el_ind = find(numbering.electrode_index == selected_el(elecID));
            rms_semggraph(el_ind)= rms_semgwindow(elecID)./1e3;        % mv --> V
         end
         rms_semggraph(rms_semggraph==0)=[];

         % rms with electrodes according to sequential order around the arm
         rms_semgreal = rms_semgwindow; 
         rms_semgRealMat = [rms_semgRealMat, rms_semgreal];

         % Estimated currents
         method = 'L2reg_area';
         ifplot = 0;
         [Vout, Iest, Vnodes] = estimateCurrents(LFM, rms_semggraph, rms_semgwindow, ...
                                        numbering.electrode_index, selected_el, graph_info.LR, method, muscleareas_norm, ifplot);

         Voutarm_f = [Vout ; Vout(1,1)];
         rms_semgreal_f = [rms_semgreal; rms_semgreal(1,1)];
         VnodesMat = [VnodesMat, Vnodes];               % voltage of the electrodes around the arm
         IestMat = [IestMat, Iest];                     % estimated currents 
         VoutMat = [VoutMat, Vout];                     % voltage of the electrodes with the graph order
         data_est.Iest = Iest;
         data_est.Voutarm_f = Voutarm_f;
         data_est.rms_semgreal_f = rms_semgreal_f;

         
         %Update of the figure
         if ifshowplot == true
             update_plot(figure_info, emg_data.time, window_start(i), semg_selected, win_time, data_est )
         end
end

%% 3)                       GOODNESS OF FIT
% =========================================================================
% Given the ouput value of the voltage estimated, assess the fitting error 
% on the real data
% =========================================================================

fprintf("The GOF for the L2-Regularized Minimum Norm Solution.\n\n")
for i = 1 : n_windows
    NRMSE_L2RArea(i) = myNRMSE(VoutMat(:,i), rms_semgRealMat(1:nelectrode,i));
end
meanNRMSE = mean(NRMSE_L2RArea);
stdNRMSE = std(NRMSE_L2RArea);

fprintf("L2 REGULARIZED MINIMUM NORM- using AREA ----------- \nGOF : %s +- %s\n\n", meanNRMSE, stdNRMSE);
        
