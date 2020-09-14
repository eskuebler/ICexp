function at = processICsweepsParFor
%{
processICsweepsParFor
- analysis of intracellular hyperpolarizing and depolarizing sweeps
%}
mainFolder = 'D:\my\';                                                      % main folder (EDIT HERE)
folder = [mainFolder,'genpath\genpath\'];                                   % general path
cellList = dir([folder,'*.mat']);                                           % list of cell data files
tic;       %10:11                                                                 % initialize clock
parfor n = 1:length(cellList)                                               % for all cells in directory
    params = loadParams;                                                    % load parameters to workspace
    cellID = cellList(n).name(1:length(cellList(n).name)-4);                % cell ID (used for saving data)
    disp(cellID)                                                            % display ID number
    a = loadFile(folder,cellList(n).name);                                  % load voltage data
    if a.LP.fullStruct == 1                                                 % if all data is present for LP
        for k = 1:length(a.LP.V)                                            % for each LP sweep
            qc = estimateRMSnoisePFfunction(a.LP,k,params,cellID,folder);   % RMS noise measurements (for QC)
            if sum(qc.logicVec) == 0                                        % if sweep passes QC criteria
                if a.LP.sweepAmps(k,1) > 0                                  % if current input > 0
                    a.LP.stats{k,1} = processDepolarizingLongPulsePF...
                        (a.LP,params,k,cellID,folder);                      % depolarizing stimulus analysis
                    a.LP.stats{k,1}.qc = qc;                                % add RMS values to data structure
                    % plateau potentials?
                else                                                        % if current input < 0
                    a.LP.stats{k,1} = processHyperpolarizingLongPulsePF...
                        (a.LP,params,qc,k,cellID,folder);                   % hyperpolarizing stimulus analysis
                    a.LP.stats{k,1}.qc = qc;                                % add RMS values to data structure
                    % variability in the steady state?
                end                                                         % end current level if
            else                                                            % if QC fails
                plotQCfailed(a.LP,k,cellID,qc,folder,params)                % plot raw voltage trace
                a.LP.stats{k,1}.qc = qc;                                    % store QC parameters
            end                                                             % end QC logic if
        end                                                                 % end sweep for loop
        LP = betweenSweepVmQC(a.LP,cellID,folder,params); a.LP = LP;        % QC between sweeps
    end                                                                     % end full structure if
    a.LP.subSummary = summarizeSubthreshold(a.LP,cellID,folder,params);     % subthreshold summary
    saveFile(a,cellID,folder);                                              % save and count cell
end                                                                         % end cell level for loop
at = toc/60;                                                                % analysis duration in seconds