function supraStats = processDepolarizingLongPulsePF(LP,params,k, ...
    cellID,folder)
%{
processSuprathresholdLongPulsePF
- analysis of depolarizing sweeps
%}
LP.putSpTimes = LP.stimOn(1,k)+...
    find(LP.V{1,k}(LP.stimOn(1,k):LP.stimOff(1,k))>=params.thresholdV)-1;   % voltage threshold
LP = getSPdVdt(LP,k,params.thresholdDVDT,cellID,folder,params);             % derivative threshold
% assess agreement betweeen detection, assign peak based on dV/dt, remove
% setting of interval for peak detection
if ~isempty(LP.putSpTimes)                                                  % if no spikes
    [int4Peak,LP.putSpTimes2] = int4APs(LP.putSpTimes);                     % interval for peak voltage
    sp = estimatePeak(LP,int4Peak,k);                                       % estimate of peak
    [sp,LP] = estimateMaxdVdtNthreshold(LP,sp,k,params,cellID,folder);      % dV/dt & initial threshold
    if ~isempty(sp.peak)                                                    % if no spikes
        [sp,LP] = refineThreshold(LP,sp,k,params);                          % refine threshold estimate
        [sp,LP] = estimateTrough(LP,sp,k,params);                           % estimate trough
        if ~isempty(sp.peak)                                                % if no spikes
            [sp,LP] = estimateSpParams(LP,sp,k,params);                     % estimate spike parameters
            if ~isempty(sp.peak)                                            % if no spikes
                sp = estimateAPTrainParams(LP,sp,k);                        % estimate spike train parameters
                % estimate plateau potential
                % assessPersistence (spikes post-stim)
                wf = getWaveforms(LP,params,sp,k);                          % get spike waveforms 
                supraStats = storeSPparams(LP,sp,wf,k);                     % store spike parameters
                plotQCdDepolarizing(LP,sp,k,cellID,folder,params)           % plot voltage and spike parameters
            else                                                            % if there are no spikes
                supraStats = outputNaNs(LP,k);                              % output structure of NaNs
            end
        else                                                                % if there are no spikes
            supraStats = outputNaNs(LP,k);                                  % output structure of NaNs
        end
    else                                                                    % if there are no spikes
        supraStats = outputNaNs(LP,k);                                      % output structure of NaNs
    end
else                                                                        % if there are no spikes
    supraStats = outputNaNs(LP,k);                                          % output structure of NaNs
end