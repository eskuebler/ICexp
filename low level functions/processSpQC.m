%{
processSpQC
%}

if isfield(a.LP.stats{k,1},'qcRemovals') && ...
        sum(isnan([ ...
        a.LP.stats{k,1}.qcRemovals.minInterval, ...
        a.LP.stats{k,1}.qcRemovals.dVdt0, ...
        a.LP.stats{k,1}.qcRemovals.mindVdt, ...
        a.LP.stats{k,1}.qcRemovals.maxThreshold, ...
        a.LP.stats{k,1}.qcRemovals.diffthreshold2peak, ...
        a.LP.stats{k,1}.qcRemovals.diffthreshold2peakT, ...
        a.LP.stats{k,1}.qcRemovals.minIntervalRe, ...
        a.LP.stats{k,1}.qcRemovals.dVdt0Re, ...
        a.LP.stats{k,1}.qcRemovals.diffpeak2trough, ...
        a.LP.stats{k,1}.qcRemovals.minTrough, ...
        a.LP.stats{k,1}.qcRemovals.percentRheobaseHeight])) ~= 11 && ...
        sum([sum(a.LP.stats{k,1}.qcRemovals.QCmatT2P), ...
        sum(a.LP.stats{k,1}.qcRemovals.QCmatT2PRe), ...
        sum(a.LP.stats{k,1}.qcRemovals.QCmatTrough), ...
        sum(a.LP.stats{k,1}.qcRemovals.percentRheobaseHeight)]) > 0
    
    vec = []; vectag = [];
    if isempty(a.LP.stats{k,1}.qcRemovals.percentRheobaseHeight)
        a.LP.stats{k,1}.qcRemovals.percentRheobaseHeight = 0;
    end
    
    figure('Position',[50 50 300 250]); set(gcf,'color','w');
%     if size(a.LP.stats{k,1}.qcRemovals.QCmatT2P,1) == 1
%         imagesc([a.LP.stats{k,1}.qcRemovals.QCmatT2P, ...
%             a.LP.stats{k,1}.qcRemovals.QCmatT2PRe, ...
%             a.LP.stats{k,1}.qcRemovals.QCmatTrough, ...
%             a.LP.stats{k,1}.qcRemovals.percentRheobaseHeight])
%     else
        a_orig = size(a.LP.stats{k,1}.qcRemovals.QCmatT2P,1);
        if sum(sum(a.LP.stats{k,1}.qcRemovals.QCmatT2PRe))==0
            a.LP.stats{k,1}.qcRemovals.QCmatT2PRe = ...
                zeros(a_orig, ...
                size(a.LP.stats{k,1}.qcRemovals.QCmatT2PRe,2));             % setting this vector to zeros (frustratingly necessary
        else
            [b,c] = size(a.LP.stats{k,1}.qcRemovals.QCmatT2PRe);
            if a_orig - b > 0
                a.LP.stats{k,1}.qcRemovals.QCmatT2PRe(b + 1: ...
                    a_orig,1:c) = zeros(a_orig - b,c);
            end
        end
        if sum(sum(a.LP.stats{k,1}.qcRemovals.QCmatTrough))==0
            a.LP.stats{k,1}.qcRemovals.QCmatTrough = ...
                zeros(a_orig, ...
                size(a.LP.stats{k,1}.qcRemovals.QCmatTrough,2));             % setting this vector to zeros (frustratingly necessary
        else
            [b,c] = size(a.LP.stats{k,1}.qcRemovals.QCmatTrough);
            if a_orig - b > 0
                a.LP.stats{k,1}.qcRemovals.QCmatTrough(b + 1: ...
                    a_orig,1:c) = zeros(a_orig - b,c);
            end
        end
        if sum(a.LP.stats{k,1}.qcRemovals.QCmatpercentRheobaseHeight)==0
            a.LP.stats{k,1}.qcRemovals.QCmatpercentRheobaseHeight = ...
                zeros(size(a.LP.stats{k,1}.qcRemovals.QCmatT2P,1),1);       % setting this vector to zeros (frustratingly necessary
        else
            [b,c] = size(a.LP.stats{k,1}.qcRemovals.QCmatpercentRheobaseHeight);
            if a_orig - b > 0
                a.LP.stats{k,1}.qcRemovals.QCmatpercentRheobaseHeight(1:c, ...
                    b + 1: a_orig) = zeros(c,a_orig - b);
            end
        end
        imagesc([a.LP.stats{k,1}.qcRemovals.QCmatT2P, ...
            a.LP.stats{k,1}.qcRemovals.QCmatT2PRe, ...
            a.LP.stats{k,1}.qcRemovals.QCmatTrough, ...
            a.LP.stats{k,1}.qcRemovals.QCmatpercentRheobaseHeight])
%     end
    colormap('gray')
    xticks(1:11)
    xticklabels({'interval','null dV/dt','dV/dt<5mV/ms', ...
        'threshold>-20mV','t2p<20mV','t2pT>2ms','interval Re', ...
        'null dV/dt Re','p2t<30mV','trough>-30mV','<30% Rheobase height'})
    xtickangle(45)
    ylabel('spike #')
    title([int2str(a.LP.sweepAmps(k,1)),' pA'])
    export_fig(['D:\my\genpath\',cellID,' ',int2str(k),' spike QC matrix'],'-pdf','-r100');
    close

    %{
    ADD other parameters below, n=11 in total all spike-wise criteria
    %}
    
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