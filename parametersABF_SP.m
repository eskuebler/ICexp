% parametersABF_SP

SP.acquireRes = h.fADCSampleInterval/1000;                                     % resolution of acquisition
if SP.acquireRes == 0.1
    constantShift = 0;                                                   % shift for input versus response
elseif SP.acquireRes == 0.05
    constantShift = 626;
    if sum(fileList(k).name(1:4)=='2018')==4 ...
            || sum(fileList(k).name(1:2)=='18')==2 ...
            || sum(fileList(k).name(1:4)=='2019')==4 ...
            || sum(fileList(k).name(1:2)=='19')==2 
        constantShift = 32;
    end
end
SP.stimOn = h.DACEpoch.lEpochInitDuration(1) + constantShift;                  % stimulus turns on
SP.stimOff = SP.stimOn+h.DACEpoch.lEpochInitDuration(2)+(postSP/SP.acquireRes);                        % stimulus turns off
for swp = 1:size(d,2)
    SP.V{1,swp} = d(SP.stimOn-(preSP/SP.acquireRes):SP.stimOff,swp)';
    SPcount = SPcount + 1;
end
SP.stimOn = SP.stimOn-(SP.stimOn-(preSP/SP.acquireRes));
SP.stimOff = SP.stimOff-(SP.stimOn-(preSP/SP.acquireRes));
clear d
tempC = length(SP.V);
SP.input = h.DACEpoch.fEpochInitLevel(2);
SP.inputInc = h.DACEpoch.fEpochLevelInc(2);
SP.sweepAmps = SP.input+(0:tempC-1)*SP.inputInc;
clear h constantShift tempC