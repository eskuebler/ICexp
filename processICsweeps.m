%{
processICsweeps
%}

clear; close all; clc;

params = loadParams;
cellList = dir('D:\genpath\genpath\*.mat');

for n = 1:length(cellList)                                                  % for all cells in directory
    tic
    cellID = cellList(n,1).name(1:length(cellList(n,1).name)-4);            % temp cell ID
    disp(cellID)
    load(['D:\genpath\genpath\',cellList(n).name])                          % load .mat files
    LP.qcSweeps = 1; LP.qcSubSweeps = 1; SP.qcSweeps = 1;                   % initialize QC'd sweep count
    if LP.fullStruct == 1                                                   % if all data is present for LP
        for k = 1:length(LP.V)                                              % for each LP sweep
            estimateRMSnoise                                                % RMS noise measurements
            if rmse_pre_st < params.RMSEst && ...
                rmse_post_st < params.RMSEst &&...
                rmse_pre < params.RMSElt && ...
                rmse_post < params.RMSElt && ...
                diffV_b_e < params.maxDiffBwBeginEnd && ...
                LP.restVPre(1,k) < params.minimumRestingPot                 % within-sweep QC
                if LP.sweepAmps(k,1) > 0                                % if suprathreshold
                    supraStats = processSuprathresholdLongPulse(LP,cellID,params,k);
                else
                    processSubthresholdLongPulse
                end
            end
            clear rmse_pre_st rmse_post_st rmse_pre rmse_post diffV_b_e diffBwSweeps
        end
    end
    LP.qcSweeps = LP.qcSweeps-1; LP.qcSubSweeps = LP.qcSubSweeps-1;
%     if SP.fullStruct == 1
%         for k = 1:length(LP.V)
%             estimateRMSnoise                                                % RMS measurements
%             if rms_noise_st < 0.07 && ...
%                     rms_noise_lt < 0.5 && ...
%                     diff_b_e < 1
%                 processShortPulse
%             end
%             clear rms_noise_st rms_noise_lt diff_b_e
%         end
%     end

    summarySubthresholdLP

    if exist('LP') && exist('SP') && exist('trees')
        save(['D:\genpath\',cellID,'.mat'],'LP','SP','trees'); clear LP SP trees
    elseif ~exist('LP') && exist('SP') && exist('trees')
        save(['D:\genpath\',cellID,'.mat'],'SP','trees'); clear SP trees
    elseif exist('LP') && ~exist('SP') && exist('trees')
        save(['D:\genpath\',cellID,'.mat'],'LP','trees'); clear LP trees
    elseif exist('LP') && exist('SP') && ~exist('trees')
        save(['D:\genpath\',cellID,'.mat'],'LP','SP'); clear LP SP
    elseif exist('LP') && ~exist('SP') && ~exist('trees')
        save(['D:\genpath\',cellID,'.mat'],'LP'); clear LP
    elseif exist('SP') && ~exist('LP') && ~exist('trees')
        save(['D:\genpath\',cellID,'.mat'],'SP'); clear SP
    end
    
    analysisT(n) = toc;
    disp(['projected time: ',int2str((round(mean(analysisT))*(length(cellList)-n)/60)/60),' hours'])
    
end