function supraStats = processSuprathresholdLongPulsePF(LP,params,k,cellID)

%{
processSuprathresholdLongPulsePF
- analysis of suprathreshold sweeps
%}

LP.putSpTimes = LP.stimOn(1,k)+...
    find(LP.V{1,k}(LP.stimOn(1,k):LP.stimOff(1,k))>=params.thresholdV)-1;   % voltage threshold
LP = getSPdVdt(LP,k,params.thresholdDVDT,cellID);                           % derivative threshold
if ~isempty(LP.putSpTimes)                                                  % if no spikes
    [int4Peak,LP.putSpTimes2] = int4APs(LP.putSpTimes);                     % interval for peak voltage
    if ~isempty(LP.putSpTimes2)                                             % if no spikes
        [sp] = estimatePeak(LP,int4Peak,k);                                 % estimate of peak
        [sp,LP] = estimateMaxdVdtNthreshold(LP,sp,k,params,cellID);         % dV/dt & initial threshold
        if ~isempty(sp.peak)                                                % if no spikes
            [sp] = refineThreshold(LP,sp,k,params);                         % refine threshold estimate
            [sp,LP] = estimateTrough(LP,sp,k,params);                       % estimate trough
            if ~isempty(sp.peak)                                            % if no spikes
                [sp,LP] = estimateSpParams(LP,sp,k,params);                 % estimate spike parameters
                if ~isempty(sp.peak)                                        % if no spikes
                    [sp] = estimateAPTrainParams(LP,sp,k);                  % estimate spike train parameters
                    % estimate plateau potential
                    % assessPersistence
                    wf = getWaveforms(LP,params,sp,k);                      % get spike waveforms 
                    supraStats = storeSPparams(LP,sp,wf,k);                 % store spike parameters
                    plotSuprathreshold(LP,sp,k,cellID)                      % plot voltage and spike parameters
                else                                                        % if there are no spikes
                    supraStats = outputNaNs(LP,k);                          % output structure of NaNs
                end
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