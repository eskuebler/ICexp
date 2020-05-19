%{
processSpQC
%}

if isfield(a.LP.stats{k,1},'qcRemovals')
    
    vec = []; vectag = [];
    LP.qcRemovals.QCmatT2P
    LP.qcRemovals.QCmatT2PRe
    LP.qcRemovals.QCmatTrough
    LP.qcRemovals.percentRheobaseHeight
    
    figure('Position',[50 50 200 250]); set(gcf,'color','w');
    imagesc(LP.qcRemovals.QCmatT2P)
    colormap('gray')
    colorbar
    xticks(1:6)
    xticklabels({'interval','null dV/dt','dV/dt<5mV/ms','threshold>-20mV','t2p<20mV','t2pT>2ms'})
    xtickangle(45)
    ylabel('spike #')
    export_fig(['D:\genpath\',cellID,' ',int2str(k),' suprathreshold spike QC'],'-pdf','-r100');
    close

    % max threshold
    if isfield(a.LP.stats{k,1}.qcRemovals,'maxThreshold') & ...
            ~isnan(a.LP.stats{k,1}.qcRemovals.maxThreshold)
        spqcmat(n,k) = spqcmat(n,k) + ...
            length(a.LP.stats{k,1}.qcRemovals.maxThreshold);
        spqcmatn(k,1) = length(a.LP.stats{k,1}.qcRemovals.maxThreshold);
        vec = a.LP.stats{k,1}.qcRemovals.maxThreshold;
        vectag = ones(1,length(a.LP.stats{k,1}.qcRemovals.maxThreshold));
    end
    % diff b/w threshold and peak
    if isfield(a.LP.stats{k,1}.qcRemovals,'diffthreshold2peak') & ...
            ~isnan(a.LP.stats{k,1}.qcRemovals.diffthreshold2peak)
        spqcmat(n,k) = spqcmat(n,k) + ...
            length(a.LP.stats{k,1}.qcRemovals.diffthreshold2peak);
        spqcmatn(k,2) = length(a.LP.stats{k,1}.qcRemovals.diffthreshold2peak);
        vec = [vec, a.LP.stats{k,1}.qcRemovals.diffthreshold2peak];
        vectag = [vectag,ones(1,length(a.LP.stats{k,1}.qcRemovals.diffthreshold2peak))+1];
    end
    % diff b/w threshold and peak T
    if isfield(a.LP.stats{k,1}.qcRemovals,'diffthreshold2peakT') & ...
            ~isnan(a.LP.stats{k,1}.qcRemovals.diffthreshold2peakT)
        spqcmat(n,k) = spqcmat(n,k) + ...
            length(a.LP.stats{k,1}.qcRemovals.diffthreshold2peakT);
        spqcmatn(k,3) = length(a.LP.stats{k,1}.qcRemovals.diffthreshold2peakT);
        vec = [vec, a.LP.stats{k,1}.qcRemovals.diffthreshold2peakT];
        vectag = [vectag,ones(1,length(a.LP.stats{k,1}.qcRemovals.diffthreshold2peakT))+2];
    end
    % diff b/w peak and trough
    if isfield(a.LP.stats{k,1}.qcRemovals,'diffpeak2trough') & ...
            ~isnan(a.LP.stats{k,1}.qcRemovals.diffpeak2trough)
        spqcmat(n,k) = spqcmat(n,k) + ...
            length(a.LP.stats{k,1}.qcRemovals.diffpeak2trough);
        spqcmatn(k,4) = length(a.LP.stats{k,1}.qcRemovals.diffpeak2trough);
        vec = [vec, a.LP.stats{k,1}.qcRemovals.diffpeak2trough];
        vectag = [vectag,ones(1,length(a.LP.stats{k,1}.qcRemovals.diffpeak2trough))+3];
    end
    % min Trough
    if isfield(a.LP.stats{k,1}.qcRemovals,'minTrough') & ...
            ~isnan(a.LP.stats{k,1}.qcRemovals.minTrough)
        spqcmat(n,k) = spqcmat(n,k) + ...
            length(a.LP.stats{k,1}.qcRemovals.minTrough);
        spqcmatn(k,5) = length(a.LP.stats{k,1}.qcRemovals.minTrough);
        vec = [vec, a.LP.stats{k,1}.qcRemovals.minTrough];
        vectag = [vectag,ones(1,length(a.LP.stats{k,1}.qcRemovals.minTrough))+4];
    end
    % percent of Rheobase Height
    if isfield(a.LP.stats{k,1}.qcRemovals,'percentRheobaseHeight') & ...
            ~isnan(a.LP.stats{k,1}.qcRemovals.percentRheobaseHeight)
        spqcmat(n,k) = spqcmat(n,k) + ...
            length(a.LP.stats{k,1}.qcRemovals.percentRheobaseHeight);
        spqcmatn(k,6) = length(a.LP.stats{k,1}.qcRemovals.percentRheobaseHeight);
        vec = [vec, a.LP.stats{k,1}.qcRemovals.percentRheobaseHeight];
        vectag = [vectag,ones(1,length(a.LP.stats{k,1}.qcRemovals.percentRheobaseHeight))+5];
    end
    
    % good before bad, bad before good
    if ~isempty(vec)
        vec = sort(vec);
            if isfield(a.LP.stats{k,1},'spTimes') && ...
                    sum(~isnan(a.LP.stats{k,1}.spTimes))>0
                spvec = abs(a.LP.stimOn(k)-a.LP.stats{k,1}.spTimes)*a.LP.acquireRes;
                vec = abs(a.LP.stimOn(k)-vec)*a.LP.acquireRes;
                sp = [spvec,vec];
                spid = [ones(1,length(spvec)),zeros(1,length(vec))];
                [~,I] = sort(sp);
                spid = spid(I);
                spqcmatnbinary(k,1:length(sp)) = sp;
                spqcmatnbinaryid(k,1:length(spid)) = spid;
                if length(unique(vectag))>1
                    'yo'
                end
                spqcvectag(k,1:length(vectag)) = vectag;
                input_current_spqc(k,1) = round(double(a.LP.sweepAmps(k,1)));
                clear sp spid B I
            end
    else            % perhaps there are good spikes?
        if isfield(a.LP.stats{k,1},'spTimes') && ...
            sum(~isnan(a.LP.stats{k,1}.spTimes))>0
            spvec = abs(a.LP.stimOn(k)-a.LP.stats{k,1}.spTimes)*a.LP.acquireRes;
            sp = spvec;
            spid = ones(1,length(spvec));
            [~,I] = sort(sp);
            spid = spid(I);
            spqcmatnbinary(k,1:length(sp)) = sp;
            spqcmatnbinaryid(k,1:length(spid)) = spid;
            input_current_spqc(k,1) = round(double(a.LP.sweepAmps(k,1)));
            spqcmatnbinary(k,length(sp)+1:end) = NaN;
            spqcmatnbinaryid(k,length(spid)+1:end) = NaN;
            spqcvectag(k,length(vectag)+1:end) = NaN;
            clear sp spid B I
        else
%             input_current_spqc(k,1) = round(double(a.LP.sweepAmps(k,1)));
        end
    end
else
%     input_current_spqc(k,1) = round(double(a.LP.sweepAmps(k,1)));
end