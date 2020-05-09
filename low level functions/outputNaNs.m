function supraStats = outputNaNs(LP,k)

% store spike-wise QC
if isfield(LP,'qcRemovals')
    if isfield(LP.qcRemovals,'maxThreshold')
        supraStats.qcRemovals.maxThreshold = LP.qcRemovals.maxThreshold;
    end
    if isfield(LP.qcRemovals,'diffthreshold2peak')
        supraStats.qcRemovals.diffthreshold2peak = LP.qcRemovals.diffthreshold2peak;
    end
    if isfield(LP.qcRemovals,'diffthreshold2peakT')
        supraStats.qcRemovals.diffthreshold2peakT = LP.qcRemovals.diffthreshold2peakT;
    end
    if isfield(LP.qcRemovals,'diffpeak2trough')
        supraStats.qcRemovals.diffpeak2trough = LP.qcRemovals.diffpeak2trough;
    end
    if isfield(LP.qcRemovals,'minTrough')
        supraStats.qcRemovals.minTrough = LP.qcRemovals.minTrough;
    end
    if isfield(LP.qcRemovals,'percentRheobaseHeight')
        supraStats.qcRemovals.maxThreshold = LP.qcRemovals.percentRheobaseHeight;
    end
else
    supraStats.qcRemovals.maxThreshold = NaN;
    supraStats.qcRemovals.diffthreshold2peak = NaN;
    supraStats.qcRemovals.diffthreshold2peakT = NaN;
    supraStats.qcRemovals.diffpeak2trough = NaN;
    supraStats.qcRemovals.minTrough = NaN;
    supraStats.qcRemovals.percentRheobaseHeight = NaN;
end

% store NaNs
supraStats.spTimes = NaN;
supraStats.spWaveforms = NaN;
supraStats.peak = NaN; 
supraStats.peakTime = NaN;
supraStats.maxdVdt = NaN;
supraStats.maxdVdtTime = NaN;
supraStats.threshold = NaN;
supraStats.thresholdTime = NaN;
supraStats.thresholdRef = NaN;
supraStats.thresholdRefTime = NaN;
supraStats.trough = NaN;
supraStats.troughTime = NaN;
supraStats.peak2trough = NaN;
supraStats.heightPT = NaN;
supraStats.halfHeightTimeUpPT = NaN;
supraStats.halfHeightTimeDownPT = NaN;
supraStats.fullWidthPT = NaN;
supraStats.heightTP = NaN;
supraStats.halfHeightTimeUpTP = NaN;
supraStats.halfHeightTimeDownTP = NaN;
supraStats.fullWidthTP = NaN;
supraStats.peakUpStroke = NaN;
supraStats.peakDownStroke = NaN;
supraStats.peakStrokeRatio = NaN;
supraStats.fastTrough = NaN;
supraStats.fastTroughDur = NaN;
supraStats.slowTrough = NaN;
supraStats.slowTroughDur = NaN;
supraStats.latency = NaN;
supraStats.meanFR50 = NaN;
supraStats.meanFR100 = NaN;
supraStats.meanFR250 = NaN;
supraStats.meanFR500 = NaN;
supraStats.meanFR750 = NaN;
supraStats.meanFR1000 = NaN;
supraStats.peakAdapt = NaN;
supraStats.ISI = NaN;
supraStats.instaRate = NaN;
supraStats.meanISI = NaN;
supraStats.cvISI = NaN;
supraStats.adaptIndex = NaN;
supraStats.adaptIndex2 = NaN;
supraStats.peakAdapt2 = NaN;
supraStats.delay = NaN;
supraStats.burst = NaN;
supraStats.sweepID = k;