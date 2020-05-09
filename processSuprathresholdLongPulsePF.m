function supraStats = processSuprathresholdLongPulsePF(LP,params,k,cellID)

%{
processSuprathresholdLongPulsePF
- analysis of suprathreshold sweeps
%}

% spikes based on 0 mV threshold
LP.putSpTimes = LP.stimOn(1,k)+...
    find(LP.V{1,k}(LP.stimOn(1,k):LP.stimOff(1,k))>=params.thresholdV)-1;   % 0 mV threshold

% spikes based on dV/dt 20mV/ms threshold
LP = getSPdVdt(LP,k,params.thresholdDVDT);

if ~isempty(LP.putSpTimes)                                                  % if no spikes
    [int4Peak,LP.putSpTimes2] = int4APs(LP.putSpTimes);                     % interval for peak voltage
    if ~isempty(LP.putSpTimes2)                                             % if no spikes
        [sp] = estimatePeak(LP,int4Peak,k);                                 % estimate of peak
        [sp,LP] = estimateMaxdVdtNthreshold(LP,sp,k,params);                % dV/dt & initial threshold
        if ~isempty(sp.peak)                                                % if no spikes
            [sp] = refineThreshold(LP,sp,k);                                % refine threshold estimate
            [sp,LP] = estimateTrough(LP,sp,k,params);                       % estimate trough
            if ~isempty(sp.peak)                                            % if no spikes
                [sp,LP] = estimateSpParams(LP,sp,k,params);                 % estimate spike parameters
                if ~isempty(sp.peak)                                        % if no spikes
                    [sp] = estimateAPTrainParams(LP,sp,k);                  % estimate spike train parameters
                    % estimate plateau potential
                    % assessPersistence
                    wf = getWaveforms(LP,params,sp,k);                      % get spike waveforms 
                    supraStats = storeSPparams(LP,sp,wf,k);                 % store spike parameters
                    plotSuprathreshold(LP,sp,k,cellID)
                else
                    supraStats = outputNaNs(LP,k);                          % output structure of NaNs
                    plotNoSP(LP,k,cellID)                                   % plot raw voltage trace
                end
            else
                supraStats = outputNaNs(LP,k);                              % output structure of NaNs
                plotNoSP(LP,k,cellID)                                       % plot raw voltage trace
            end
        else
            supraStats = outputNaNs(LP,k);                                  % output structure of NaNs
            plotNoSP(LP,k,cellID)                                           % plot raw voltage trace
        end
    else
        supraStats = outputNaNs(LP,k);                                      % output structure of NaNs
        plotNoSP(LP,k,cellID)                                               % plot raw voltage trace
    end
else
    supraStats = outputNaNs(LP,k);                                          % output structure of NaNs
    plotNoSP(LP,k,cellID)                                                   % plot raw voltage trace
end