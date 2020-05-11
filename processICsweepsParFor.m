function at = processICsweepsParFor
%{
processICsweepsParFor
- analysis of intracellular subthreshold and suprathreshold sweeps
%}
folder = 'D:\genpath\genpath\'; cellList = dir([folder,'*.mat']);           % list of cell data files
tic;                                                                        % initialize clock
for n = 1:length(cellList)                                                  % for all cells in directory
    params = loadParams;                                                    % load parameters to workspace
    cellID = cellList(n).name(1:length(cellList(n).name)-4);                % cell ID (used for saving data)
    disp(cellID)                                                            % display ID number
    a = loadFile(folder,cellList(n).name);                                  % load voltage data
    if a.LP.fullStruct == 1                                                 % if all data is present for LP
        for k = 1:length(a.LP.V)                                            % for each LP sweep
            qc = estimateRMSnoisePFfunction(a.LP,k,params,cellID);          % RMS noise measurements (for QC)
            if sum(qc.logicVec) == 0                                        % if QC passes
                if a.LP.sweepAmps(k,1) > 0                                  % if current input > 0
                    a.LP.stats{k,1} = processSuprathresholdLongPulsePF...
                        (a.LP,params,k,cellID); a.LP.stats{k,1}.qc = qc;    % suprathreshold analysis
                    % plateau potentials?
                else                                                        % if current input < 0
                    a.LP.stats{k,1} = processSubthresholdLongPulsePF...
                        (a.LP,params,qc,k,cellID); a.LP.stats{k,1}.qc = qc; % subthreshold analysis
                    % variability in the steady state?
                end                                                         % end current level if
            else                                                            % if QC fails
                plotQCfailed(a.LP,k,cellID,qc)                              % plot raw voltage trace
                a.LP.stats{k,1}.qc = qc;                                    % store QC parameters
            end                                                             % end QC logic if
        end                                                                 % end sweep for loop
        LP = betweenSweepVmQC(a.LP,cellID); a.LP = LP;                      % QC between sweeps
    end                                                                     % end full structure if
    a.LP.subSummary = summarizeSubthreshold(a.LP,cellID);                   % subthreshold summary
    saveFile(a,cellID);                                                     % save and count cell
end                                                                         % end cell level for loop
at = toc/60;                                                                % analysis duration in seconds