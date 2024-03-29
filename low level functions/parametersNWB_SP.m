% parametersNWB_SP

temp = h5read(fileName,[level.Resp,'/data'])'*1000;
stimFun = h5read(fileName,[level.Stim,'/data']);
if sum(h5readatt(fileName,[level.Stim,'/starting_time'],'unit') == 'Seconds')==7
    SP.acquireRes = 1000/h5readatt(fileName,[level.Stim,'/starting_time'],'rate');
end
temp2 = find(stimFun(10000:end,1)~=0, 1,'first')+10000-1;                                           % offset by 10,000 timepoints to avoid test pulses
if ~isempty(temp2)
    SP.stimOn(1,SPcount) = temp2;
    SP.stimOff(1,SPcount) = find(stimFun~=0, 1,'last')+1+(postSP/SP.acquireRes);
    stimDur = ((SP.stimOff(1,SPcount)-SP.stimOn(1,SPcount))*SP.acquireRes)-postSP;
    if length(unique(stimFun))==3 && ...
            length(findpeaks(stimFun))==2 && ...
            stimDur==3 && ...
            length(temp)>(SP.stimOff(1,SPcount)+(postSP/SP.acquireRes))                 % errors in curation
        SP.V{1,SPcount} = temp(SP.stimOn(1,SPcount)-(preSP/SP.acquireRes):SP.stimOff(1,SPcount));
        SP.stimOn(1,SPcount) = SP.stimOn(1,SPcount)-(SP.stimOn(1,SPcount)-(preSP/SP.acquireRes));
        SP.stimOff(1,SPcount) = SP.stimOff(1,SPcount)-(SP.stimOn(1,SPcount)-(preSP/SP.acquireRes));
        SP.sweepAmps(SPcount,1) = h5read(fileName,[level.Stim,'/aibs_stimulus_amplitude_pa']);
        SPcount = SPcount + 1;
    end
end
clear stimFun temp stimDur