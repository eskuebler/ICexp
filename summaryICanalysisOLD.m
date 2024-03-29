%{
summaryICanalysis
%}

clear; close all; clc;                                                      % prepare Matlab workspace

% list of cells
mainFolder = 'D:\Allen Institute Data\patchclamp\my\';                                                          % main folder for data (EDIT HERE)
cellList = dir([mainFolder,'genpath\*.mat']);                               % cell list
    
% free parameters (across sweep QC)
minNmaxThres = 10;                                                          % threshold for difference b/w min and max voltage
origStdThresMax = 3;                                                        % above this threshold we remove sweeps > 1.75 S.D.s
origStdThresMin = 0.7;                                                      % below this threshold cells are not QC'd sweep-wise any further

% summary file names
% savefilename = 'cell type details (2020 06 08)';                            % main data file name
% savefilenameInd = 'cell type details (2020 06 08) sp qc';                   % single cells file name (cell ID added below in loop)      

% load data
load('manual_entry_data.mat')                                               % load Michelle's NHP data
layerID = ID; dendrite_typeMJ = dendrite_type; clear ID dendrite_type       % adjust variable names for overlap
load('pyramidal cell IDs.mat')                                              % load Michelle's pyramidal cells
pyrID = ID; clear ID                                                        % adjust variable names for overlap
load('cell_types_specimen_details.mat','donor__species1','specimen__id1',...
    'line_name1','tag__dendrite_type1','ef__ri1','ef__tau1',...
    'ef__threshold_i_long_square1','cell_reporter_status1',...
    'structure__acronym1','structure__layer1','csl__normalized_depth');     % load Allen Institute parameters

% initialize
IC.input_current_s = zeros(length(cellList),75);
IC.firing_rate_s = zeros(length(cellList),75);
removedListStd = []; rmvdStdCount = 1;
removedListMinMax = []; rmvdMMCount = 1;
spqcmat = zeros(length(cellList),60);

for n = 1:length(cellList)                                                  % for each cells
    clc; disp(n)                                                            % display n value
    sweepIDcount = 1;
    
    % Determine from which institution the cell was obtained and assigning their matadata and control parameters
    if length(cellList(n).name) == 17
        cellID = cellList(n).name(1:length(cellList(n).name)-4);            % JMT lab ID
    else
        cellID = cellList(n).name(10:length(cellList(n).name)-4);           % Allen ID
    end

    
    % pull data for AIBS cell
    if length(cellID) == 9                                                  % if an Allen cell
        % match species ID # to get Allen parameters
        k  = find(specimen__id1==str2num(cellID));

        % get Allen Institute parameters
        IC.ID(n,:) = ([num2str(specimen__id1(k)),'__AI']);
        IC.specimen(n,1) = donor__species1(k);
        IC.struct(n,1) = structure__acronym1(k);
        IC.cortical_layer(n,1) = categorical(structure__layer1(k));
        IC.cort_depth(n,1) = csl__normalized_depth(k);
        IC.transline(n,1) = line_name1(k);
        IC.reporterStatus(n,1) = cell_reporter_status1(k);
        IC.dendrite_type(n,1) = tag__dendrite_type1(k);
        IC.AIr(n,1) = ef__ri1(k);
        IC.AItau(n,1) = ef__tau1(k);
        IC.AIrheobase(n,1) = round(ef__threshold_i_long_square1(k));
        if sum(ismember(pyrID,str2num(IC.ID(n,1:9))))
            IC.pyramidalID(n,1) = 1;
        else
            IC.pyramidalID(n,1) = 0;
        end
        IC.access_resistance(n,1) = NaN;
        IC.temperature(n,1) = 0;
        
    % pull data for PCTD cell
    else                                                                    % if a JMT cell
        IC.ID(n,:) = cellID;
        IC.specimen(n,:) = categorical(cellstr('NHP'));
        IC.struct(n,:) = categorical(cellstr('PFC'));
        IC.transline(n,:) = categorical(cellstr('N/A'));
        IC.reporterStatus(n,:) = categorical(cellstr('N/A'));
        IC.AIr(n,1) = NaN;
        IC.AItau(n,1) = NaN;
        IC.AIrheobase(n,1) = NaN;
        IC.pyramidalID(n,1) = NaN;
        if sum(ismember(layerID,cellID))
            k = find(ismember(layerID,cellID)==1);
%             if AccesResistance(k,1) ~= 'N/A'                              
%                 IC.access_resistance(n,1) = ...
%                     str2num(char(AccesResistance(k,1)));
%             end
            IC.dendrite_type(n,1) = dendrite_typeMJ(k,1);
            IC.cortical_layer(n,1) = Layer(k,1);
%             if Temperature(k,1) ~= 'N/A'
%                 IC.temperature(n,1) = str2num(char(Temperature(k,1)));
%             end
        else
            IC.access_resistance(n,1) = NaN;
            IC.dendrite_type(n,1) = categorical(cellstr('N/A'));
            IC.cortical_layer(n,1) = categorical(cellstr('N/A'));
            IC.temperature(n,1) = 0;
        end
    end

    % load our analysis
    load([mainFolder,'genpath\',cellList(n).name]);                         % load analysis file for cell
    
    % qc stuff
    qc_logic = zeros(1,6);                                                  % initialize QC matrix

    if a.LP.fullStruct == 1                                                 % if full data structure is available
        qcID{n,1} = cellID;                                                 % cell ID
        qc_V_vec(n,1:length(a.LP.rmp(1,:))) = round(a.LP.rmp(1,:),2);       % resting membrane potential
        qc_V_vecDelta(n,1:length(a.LP.rmp(1,:))) = ...
            round(a.LP.rmp(1,1) - a.LP.rmp(1,:),2);                         % diff RMP pre and post stimulus
        spqcmatn = zeros(length(a.LP.sweepAmps),10);                        % initialize count of QC removals matrix (each column is a criteria)
        binaryMatCount = 1;
        spqcvectag = nan(20,300);                                           % initialize QC tag storage
        input_current_spqc = zeros(20,1);                                   % initialize input current storage
        for k = 1:length(a.LP.sweepAmps)                                    % for each sweep
            
            % sweep-wise quality control parameters
            qc_restVpre(n,k) = round(a.LP.stats{k,1}.qc.restVPre,2);        % RMP pre stimulus
            qc_restVpost(n,k) = round(a.LP.stats{k,1}.qc.restVPost,2);      % RMP post stimulus
            qc_restVdiffpreNpost(n,k) = round( ...
                a.LP.stats{k,1}.qc.diffV_b_e,2);                            % diff RMP pre and post stimulus
            qc_rmse_pre_lt(n,k) = round(a.LP.stats{k,1}.qc.rmse_pre,2);     % long term RMS pre stimulus
            qc_rmse_post_lt(n,k) = round(a.LP.stats{k,1}.qc.rmse_post,2);   % lt RMS post stimulus
            qc_rmse_pre_st(n,k) = round(a.LP.stats{k,1}.qc.rmse_pre_st,2);  % st RMS pre stimulus
            qc_rmse_post_st(n,k) = round( ...
                a.LP.stats{k,1}.qc.rmse_post_st,2);                         % st RMS post stimulus
            qc_logic = qc_logic+a.LP.stats{k,1}.qc.logicVec;                % QC logic vector (each column is a criteria)
            
            % spike-wise QC processing            
            processSpQCOLD                                                     % process spike-wise QC
            % remove sweeps that exceed good/bad spike ratio  0.3
            % add to QC vec
            
            % assess the removal of this sweep
            if sum(a.LP.stats{k,1}.qc.logicVec) == 0 %isfield(a.LP,'stats') &&k<=length(a.LP.stats) && ...~isempty(a.LP.stats{k,1}) && ...
                qc.sweepID(n,k) = k;
                sweepBinary(n,k) = 1;
                sweepBinaryOrig(1,k) = 1;
            else
                qc.sweepID(n,k) = 0;
                sweepBinary(n,k) = 0;
                sweepBinaryOrig(1,k) = 0;
            end
        end
        if sum(sum(spqcmatn))>0                                             % figure of spike-wise QC (sweeps X criteria)
            figure('Position',[50 50 300 300]); set(gcf,'color','w');
            imagesc(spqcmatn)
            box off
            colorbar
            colormap('gray')
            xticks(1:10)
            xticklabels({'interval','null dV/dt','dV/dt<5mV/ms', ...
                'threshold>-20mV','t2pN(B)<35(45)mV','t2pT>1.5ms','interval Re', ...
                'null dV/dt Re','trough>-30mV','<30% Rheobase height'})
            xtickangle(45)
            ylabel('current input (pA)')
            yticks(1:size(spqcmatn,1))
            yticklabels({a.LP.sweepAmps})
%             export_fig(['D:\my\genpath\',cellID,' spike QC (sweeps by critiera)'],'-pdf','-r100');
            close
        end 
        if exist('spqcmatnbinaryid','var')                                  % figure of spike-wise QC (sweeps X spikes)
            figure('Position',[50 50 300 250]); set(gcf,'color','w');
            imagesc(spqcmatnbinaryid)
            box off
            colormap('gray')
            xlabel('spike # (white==passes QC)')
            ylabel('current input (pA)')
            yticks(1:length(k_len_spID))
            yticklabels({k_len_spID})
%             export_fig(['D:\my\genpath\',cellID,' spike QC (sweeps by spike binary)'],'-pdf','-r100');
            close
            clear spqcmatnbinaryid spqcmatnbinary k_len_spID
        end
        % add sp qc
        qc_logic_mat(n,1:6) = qc_logic;
        processBwSweepsQCOLD                                                   % across sweep QC
    end
    if size(qc_class_mat,1)~=n
        qc_class_mat(n,:) = 0;
    end
    
    if size(sweepBinary,1)==n && sum(sweepBinary(n,:))>0
        
        % subthreshold summary parameters
        IC.acquireRes(n,1) = double(a.LP.acquireRes);
        IC.resistance(n,1) = round(double(a.LP.subSummary.resistance),2);   % Allen Method
        IC.resistance_hd(n,1) = resistance_hd(a.LP);                        % resistance (HD?)
        IC.resistance_ss(n,1) = resistance_ss(a.LP);                        % resistance based on steady state
        IC.Vrest(n,1) = restingMP(a.LP);                                    % resting membrane potential
        IC.time_constant(n,1) = round(double(a.LP.subSummary.tauMin),2);
        IC.time_constant2(n,1) = round(double(a.LP.subSummary.tauSS),2);
        if a.LP.fullStruct == 1
            k = find(ismember(qc.sweepID(n,:),...
                find(round(double(a.LP.sweepAmps)) == -90))==1);            % find -90 pA input
            if length(k)>1  
                k = k(1);
            end
            if ~isempty(k)
                getSubthresholdStats                                        % get subthreshold stats
            else                                                            % if no -90 pA sweep
                k = find(ismember(qc.sweepID(n,:),...
                    find(round(double(a.LP.sweepAmps)) == -70))==1);        % find -70 pA sweep
                if length(k)>1
                    k = k(1);
                end
                if ~isempty(k)
                    getSubthresholdStats                                    % get subthreshold stats
                else                                                        % if no -70 pA sweep
                    k = find(ismember(qc.sweepID(n,:),...
                        find(round(double(a.LP.sweepAmps)) == -110))==1);   % find -110 pA sweep
                    if length(k)>1
                        k = k(1);
                    end
                    if ~isempty(k)
                        getSubthresholdStats                                % get subthreshold stats
                    else                                                    % if no -50 pA sweeps
                        IC.subamp(n,1) = NaN;                               % add NaNs (blank spaces in csv format)
                        IC.submin(n,1) = NaN;
                        IC.rebound_slope(n,1) = NaN;
                        IC.rebound_depolarization(n,1) = NaN;
                        IC.nb_rebound_sp(n,1) = 0;
                        IC.sag(n,1) = NaN;
                        IC.steadystate(n,1) = NaN;
                        IC.sag_ratio(n,1) = NaN;
                    end
                end
            end

            % find rheobase and parameters of first spike
            [B,I] = sort(round(double(a.LP.sweepAmps(1:length(a.LP.stats)))));
            int_vec = find(B>0);
            temp = int_vec(find(ismember(int_vec,qc.sweepID(n,:))==1));
            idx = I(temp);
            amp = B(temp);
            IC.input_current_s(n,1:length(amp)) = round(double(a.LP.sweepAmps(idx)));
            
            spCheck = 0;
            for k = 1:length(idx)
                if isfield(a.LP.stats{idx(k),1},'spTimes') && ...
                        sum(~isnan(a.LP.stats{idx(k),1}.spTimes))>0
                    % spike train parameters
                    IC.rate_1s(n,k) = a.LP.stats{idx(k),1}.meanFR1000;
                    IC.rate_750ms(n,k) = a.LP.stats{idx(k),1}.meanFR750;
                    IC.rate_500ms(n,k) = a.LP.stats{idx(k),1}.meanFR500;
                    IC.rate_250ms(n,k) = a.LP.stats{idx(k),1}.meanFR250;
                    IC.rate_100ms(n,k) = a.LP.stats{idx(k),1}.meanFR100;
                    IC.rate_50ms(n,k) = a.LP.stats{idx(k),1}.meanFR50;
                    IC.delay(n,k) = round(double(a.LP.stats{idx(k),1}.delay),2);
                    IC.burst(n,k) = round(double(a.LP.stats{idx(k),1}.burst),2);
                    IC.latency(n,k) = round(double(a.LP.stats{idx(k),1}.latency),2);
                    IC.cv_ISI(n,k) = round(double(a.LP.stats{idx(k),1}.cvISI),2);
                    IC.adaptation1(n,k) = round(double(a.LP.stats{idx(k),1}.adaptIndex),2);
                    IC.adaptation2(n,k) = round(double(a.LP.stats{idx(k),1}.adaptIndex2),2);
                    IC.peak_adaptation1(n,k) = round(double(a.LP.stats{idx(k),1}.peakAdapt),2);
                    IC.peak_adaptation2(n,k) = round(double(a.LP.stats{idx(k),1}.peakAdapt2),2);
                    spCheck = 1;
                end
            end

            for k = 1:length(idx)
                if isfield(a.LP.stats{idx(k),1},'spTimes') && ...
                        sum(~isnan(a.LP.stats{idx(k),1}.spTimes))>0
                    IC.rheobaseLP(n,1) = amp(k);
                    IC.rate_1sRheobase(n,1) = a.LP.stats{idx(k),1}.meanFR1000;
                    IC.rate_750msRheobase(n,1) = a.LP.stats{idx(k),1}.meanFR750;
                    IC.rate_500msRheobase(n,1) = a.LP.stats{idx(k),1}.meanFR500;
                    IC.rate_250msRheobase(n,1) = a.LP.stats{idx(k),1}.meanFR250;
                    IC.rate_100msRheobase(n,1) = a.LP.stats{idx(k),1}.meanFR100;
                    IC.rate_50msRheobase(n,1) = a.LP.stats{idx(k),1}.meanFR50;
                    IC.delayRheobase(n,1) = round(double(a.LP.stats{idx(k),1}.delay),2);
                    IC.burstRheobase(n,1) = round(double(a.LP.stats{idx(k),1}.burst),2);
                    IC.latencyRheobase(n,1) = round(double(a.LP.stats{idx(k),1}.latency),2);
                    IC.cv_ISIRheobase(n,1) = round(double(a.LP.stats{idx(k),1}.cvISI),2);
                    IC.adaptation1Rheobase(n,1) = round(double(a.LP.stats{idx(k),1}.adaptIndex),2);
                    IC.adaptation2Rheobase(n,1) = round(double(a.LP.stats{idx(k),1}.adaptIndex2),2);
                    IC.peak_adaptation1Rheobase(n,1) = round(double(a.LP.stats{idx(k),1}.peakAdapt),2);
                    IC.peak_adaptation2Rheobase(n,1) = round(double(a.LP.stats{idx(k),1}.peakAdapt2),2);
                    IC.peakLP(n,1) = round(double(a.LP.stats{idx(k),1}.peak(1)),2);
                    IC.peakTimeLP(n,1) = double(a.LP.stats{idx(k),1}.peakTime(1));
                    IC.thresholdLP(n,1) = round(double(a.LP.stats{idx(k),1}.thresholdRef(1)),2);
                    IC.half_width_threshold_peak(n,1) = round(double(a.LP.stats{idx(k),1}.fullWidthTP(1)),2);
                    IC.half_width_peak_trough(n,1) = round(double(a.LP.stats{idx(k),1}.fullWidthPT(1)),2);
                    IC.height_threshold_peak(n,1) = round(double(a.LP.stats{idx(k),1}.heightTP(1)),2);
                    IC.height_peak_trough(n,1) = round(double(a.LP.stats{idx(k),1}.heightPT(1)),2);
                    IC.peak_up_stroke(n,1) = round(double(a.LP.stats{idx(k),1}.peakUpStroke(1)),2);
                    IC.peak_down_stroke(n,1) = round(double(a.LP.stats{idx(k),1}.peakDownStroke(1)),2);
                    IC.peak_stroke_ratio(n,1) = round(double(a.LP.stats{idx(k),1}.peakStrokeRatio(1)),2);
                    IC.trough(n,1) = round(double(a.LP.stats{idx(k),1}.trough(1)),2);
                    IC.fastTrough(n,1) = a.LP.stats{idx(k),1}.fastTroughDur(1);
                    IC.slowTrough(n,1) = a.LP.stats{idx(k),1}.slowTroughDur(1);
                    IC.fastTroughV(n,1) = a.LP.stats{idx(k),1}.fastTrough(1);
                    IC.slowTroughV(n,1) = a.LP.stats{idx(k),1}.slowTrough(1);
                    if IC.slowTrough(n,1) < 2                               %  check if the cell has a true fAHP, its peak should be more depolarized than any other afterhyperpolarizations                         
                        IC.fAHPamp(n,1) = IC.thresholdLP(n,1) - ...
                            IC.trough(n,1);                                 %  amplitude of fAHP from threshold
                    else 
                        IC.fAHPamp(n,1) = NaN;
                    end
                    if size(a.LP.stats{idx(k),1}.waves,1) > 1
                        IC.wf(n,:) = round(mean(a.LP.stats{idx(k),1}.waves),2);
                    else
                        IC.wf(n,:) = round(a.LP.stats{idx(k),1}.waves,2);
                    end
                    disp([IC.fastTrough(n,1),IC.slowTrough(n,1)])
                    disp([a.LP.stats{idx(k),1}.fastTrough(1),a.LP.stats{idx(k),1}.slowTrough(1)])
                    figure('Position',[50 50 300 250]); set(gcf,'color','w');
                    plot(a.LP.V{1,idx(k)},'k')
                    hold on
                    scatter(a.LP.stimOn(1,idx(k))+IC.peakTimeLP(n,1),IC.peakLP(n,1),30,'g')
                    scatter(a.LP.stimOn(1,idx(k))+IC.peakTimeLP(n,1)+(IC.fastTrough(n,1)/a.LP.acquireRes),a.LP.stats{idx(k),1}.fastTrough(1),30,'r')
                    scatter(a.LP.stimOn(1,idx(k))+IC.peakTimeLP(n,1)+(IC.slowTrough(n,1)/a.LP.acquireRes),a.LP.stats{idx(k),1}.slowTrough(1),50,'k')
                    xlabel('time')
                    ylabel('voltage (mV)')
                    axis tight
                    xlim([a.LP.stimOn(1,idx(k))+IC.peakTimeLP(n,1)-(1/a.LP.acquireRes) ...
                        a.LP.stimOn(1,idx(k))+IC.peakTimeLP(n,1)+((IC.slowTrough(n,1)+1)/a.LP.acquireRes)])
                    box off
                    annotation('textbox',[.3 .7 0 0],'String',IC.fastTrough(n,1),'FitBoxToText','on','LineStyle','none','color','r');
                    annotation('textbox',[.3 .9 0 0],'String',IC.slowTrough(n,1),'FitBoxToText','on','LineStyle','none');
                    annotation('textbox',[.52 .7 0 0],'String',IC.fastTroughV(n,1),'FitBoxToText','on','LineStyle','none','color','r');
                    annotation('textbox',[.52 .9 0 0],'String',IC.slowTroughV(n,1),'FitBoxToText','on','LineStyle','none');
                    export_fig(['D:\figs\',cellID,' ',int2str(idx(k)),' trough'],'-pdf','-r300');
                    close
                    
                    break
                else
                    IC.rheobaseLP(n,1) = NaN;
                end
            end
            
            % global spike parameters and Hero sweep selection
            k = [];                                                         % resetting k for indexing sweeps
            flag = 0;                                                       % variable to fire the if condition in while loop only one time
            if spCheck == 1 
                
                % global spiketrain parameters
                IC.maxFiringRate(n,1) = max(IC.rate_1s(n,:));
                [IC.mdn_insta_freq(n,1),IC.mdn_insta_freq_times{n,1}] ...
                    = median_isi(a.LP);                                     % obtain the median ISI of all suprathreshold sweeps
                
                % picking "Hero sweep" for more spike train parameters per cell
                [~,k] = min(abs(double(B)-(IC.rheobaseLP(n,1)*1.5)));        % hero sweep is 1.5x Rheobase
                if k > 0                                                    % if there is a sweep 1.5x Rheobase
                    while ~ismember(qc.sweepID(n,k),k) ||    ...            % Making sure the k sweep meets other necessary conditions
                            ~isfield(a.LP.stats{k,1},'burst')  ||   ...     % It has spike train analysis fields like burst
                            a.LP.sweepAmps(k) <= IC.rheobaseLP(n,1) || ...  % It is not lower than the rheobase
                            a.LP.sweepAmps(k) > 3*IC.rheobaseLP(n,1)        % It is not more than triple the rheobase 
                        k = k - 1;  
                        if k == 0 && flag == 0
                            [~, k] = min(abs(double(B) - ...
                                (IC.rheobaseLP(n,1)*3)));                   % hero sweep is 8x Rheobase sweep
                            flag = 1;                                       % set if condition to fire
                        end
                        if k == 0 && flag == 1; break; end
                    end  
                end
                if k == 0 
%                     IC.rate_s_hero(n,1) = NaN;
                    IC.rate_1sHero(n,1) = NaN;
                    IC.rate_750msHero(n,1) = NaN;
                    IC.rate_500msHero(n,1) = NaN;
                    IC.rate_250msHero(n,1) = NaN;
                    IC.rate_100msHero(n,1) = NaN;
                    IC.rate_50msHero(n,1) = NaN;
                    IC.rateIntWindows(n,1:6) = NaN;
                    IC.burst_hero(n,1) = NaN;
                    IC.delay_hero(n,1) =   NaN;
                    IC.latency_hero(n,1) = NaN;
                    IC.cv_ISI(n,1) =  NaN;
                    IC.adaptation2(n,1) = NaN;
                    IC.hero_amp(n,1) = NaN;
                    IC.peakTimeLPHero{n,1} = NaN;
                elseif length(k) > 1
%                     IC.rate_s_hero(n,1) = mean(IC.rate_1s(n,find(IC.input_current_s(n,:)==a.LP.sweepAmps(k))));
                    IC.rate_1sHero(n,1) = mean(a.LP.stats{k,1}.meanFR1000);
                    IC.rate_750msHero(n,1) = mean(a.LP.stats{k,1}.meanFR750);
                    IC.rate_500msHero(n,1) = mean(a.LP.stats{k,1}.meanFR500);
                    IC.rate_250msHero(n,1) = mean(a.LP.stats{k,1}.meanFR250);
                    IC.rate_100msHero(n,1) = mean(a.LP.stats{k,1}.meanFR100);
                    IC.rate_50msHero(n,1) = mean(a.LP.stats{k,1}.meanFR50);
                    IC.rateIntWindows(n,:) = a.LP.stats{k,1}.frIntWindows;
                    IC.burst_hero(n,1) = mean(train_burst(n,k(1:length(k))));
                    IC.delay_hero(n,1) =   mean(train_delay(n,k(1:length(k))));
                    IC.latency_hero(n,1) = mean(train_latency(n,k(1:length(k))));
                    IC.cv_ISI(n,1) =   mean(train_cv_ISI(n,k(1:length(k))));
                    IC.adaptation1(n,1) = mean(train_adaptation1(n,k(1:length(k))));
                    IC.adaptation2(n,1) = mean(train_adaptation2(n,k(1:length(k))));
                    IC.hero_amp(n,1) = unique(a.LP.sweepAmps(k));     
                    IC.peakTimeLPHero{n,1} = double(a.LP.stats{k,1}.peakTime);
                else
%                     IC.rate_s_hero(n,1) = IC.rate_1s(n,find(IC.input_current_s(n,:)==a.LP.sweepAmps(k)));
                    IC.rate_1sHero(n,1) = a.LP.stats{k,1}.meanFR1000;
                    IC.rate_750msHero(n,1) = a.LP.stats{k,1}.meanFR750;
                    IC.rate_500msHero(n,1) = a.LP.stats{k,1}.meanFR500;
                    IC.rate_250msHero(n,1) = a.LP.stats{k,1}.meanFR250;
                    IC.rate_100msHero(n,1) = a.LP.stats{k,1}.meanFR100;
                    IC.rate_50msHero(n,1) = a.LP.stats{k,1}.meanFR50;
                    IC.rateIntWindows(n,:) = a.LP.stats{k,1}.frIntWindows;
                    IC.burst_hero(n,1) =  a.LP.stats{k, 1}.burst;
                    IC.delay_hero(n,1) =   a.LP.stats{k, 1}.delay;
                    IC.latency_hero(n,1) = a.LP.stats{k, 1}.latency;
                    IC.cv_ISI(n,1) =   a.LP.stats{k, 1}.cvISI ;  
                    IC.adaptation1(n,1) = a.LP.stats{k, 1}.adaptIndex; 
                    IC.adaptation2(n,1) = a.LP.stats{k, 1}.adaptIndex2; 
                    IC.hero_amp(n,1) = a.LP.sweepAmps(k);
                    IC.peakTimeLPHero{n,1} = double(a.LP.stats{k,1}.peakTime);
                end    
                clear B I idx amp temp k
            end
        end
    end
end

IC.cellList = cellList;

clear cell_reporter_status donor__species line_name specimen__id structure__acronym ...
    structure__layer tag__dendrite_type

% Cleaning up variables 
IC.resistance_ss(IC.resistance_ss==0)= NaN;                                 % 0 is a non-sensical value in most variables and should be NaN
IC.resistance_hd(IC.resistance_hd==0)= NaN;
IC.time_constant(IC.time_constant==0)= NaN;
IC.Vrest_sag_sweep(IC.Vrest_sag_sweep==0)= NaN;                           
IC.sag_ratio(IC.sag_ratio==0)= NaN;                                        
IC.mdn_insta_freq(IC.mdn_insta_freq==0)= NaN;  
IC.rheobaseLP(IC.rheobaseLP==0)= NaN;  
IC.hero_amp(IC.hero_amp==0)=NaN;
IC.firing_rate_s(IC.firing_rate_s==0)= NaN;

removeFailedQC
% ind = find(diffMinMaxV~=0);
% figure('Position',[50 50 300 250]); set(gcf,'color','w');
% histogram(diffMinMaxV(ind),40,'FaceColor','k');
% line([minNmaxThres,minNmaxThres],[0,30], ...
%             'color','r','linewidth',1,'linestyle','--');
% xlabel('diff b/w min and max V')
% ylabel('probability')
% axis tight
% box off
% % export_fig(['qc mean diff min max resting V'],'-pdf','-r100');
% % close
% 
% ind = find(meanOrigV~=0);
% figure('Position',[50 50 300 250]); set(gcf,'color','w');
% histogram(meanOrigV(ind),40,'FaceColor','k');
% xlabel('mean in V across sweep')
% ylabel('probability')
% axis tight
% box off
% % export_fig(['qc mean resting V'],'-pdf','-r100');
% % close
% 
% ind = find(stdOrigV~=0);
% figure('Position',[50 50 300 250]); set(gcf,'color','w');
% hold on
% histogram(stdOrigV(ind),40,'FaceColor','k');
% line([origStdThresMax,origStdThresMax],[0,20], ...
%             'color','r','linewidth',1,'linestyle','--');
% line([origStdThresMin,origStdThresMin],[0,20], ...
%             'color','g','linewidth',1,'linestyle','--');
% xlabel('std in V across sweep')
% ylabel('probability')
% axis tight
% box off
% % export_fig(['qc std resting V'],'-pdf','-r100');
% close

% QCanalysis
%{
% all parameters
T = table(ID,specimen,struct,cortical_layer,transline,reporterStatus,...
    dendrite_type,pyramidalID,temperature,resistance,access_resistance,...
    time_constant,time_constant2,subamp,submin,rebound_slope,nb_rebound_sp,sag,...
    steadystate,sag_ratio,rheobaseLP,delay,burst,latency,cv_ISI,adaptation1,...
    adaptation2,peak_adaptation1,peak_adaptation2,peakLP,thresholdLP,...
    half_width_threshold_peak,half_width_peak_trough,...
    height_threshold_peak,height_peak_trough,...
    peak_up_stroke,peak_down_stroke,peak_stroke_ratio,trough,...
    fastTrough,slowTrough,maxFiringRate);
    writetable(T,[savefilename,'.xlsx'],'Sheet','Sheet1',...
        'WriteRowNames',true)
% sweeps that pass qc
T = table(ID,qc_class_mat);
    writetable(T,[savefilename,'.xlsx'],'Sheet','Sheet2',...
        'WriteRowNames',true)
% spike trains all sweeps
T = table(ID,input_current_s,firing_rate_s,train_delay,train_burst,...
    train_latency,train_cv_ISI,train_adaptation1,train_adaptation2,...
    train_peak_adaptation1,train_peak_adaptation2);
    writetable(T,[savefilename,'.xlsx'],'Sheet','Sheet3',...
        'WriteRowNames',true)
% qc parameters for all sweeps
T = table(ID,meanOrigV,diffMinMaxV,stdOrigV,qc_V_vec,qc_restVpre,qc_restVpost,qc_restVdiffpreNpost,...
    qc_rmse_pre_lt,qc_rmse_post_lt,qc_rmse_pre_st,qc_rmse_post_st);
    writetable(T,[savefilename,'.xlsx'],'Sheet','Sheet4',...
        'WriteRowNames',true)
% waveforms
T = table(ID,wf);
    writetable(T,[savefilename,'.xlsx'],'Sheet','Sheet5',...
        'WriteRowNames',true)
% close excel
system('taskkill /F /IM EXCEL.EXE');
%}

% plot our analysis versus Allen Institute
% listnan = isnan(resistance);
% resistance(listnan) = [];
% AIr(listnan) = [];
% listnan = isnan(time_constant);
% time_constant(listnan) = [];
% AItau(listnan) = [];
% listnan = isnan(rheobaseLP);
% rheobaseLP(listnan) = [];
% AIrheobase(listnan) = [];

figure('Position',[50 50 450 120]); set(gcf,'color','w');
subplot(1,3,1)
scatter(IC.AIr,IC.resistance,5,'k')
xlabel('resistance (AI)')
ylabel('resistance')
[r_RI,p_Ri] = corr(IC.AIr,IC.resistance, 'rows','complete');
subplot(1,3,2)
scatter(IC.AItau,IC.time_constant,5,'k')
xlabel('time constant (AI)')
ylabel('time constant')
[r_tau,p_tau] = corr(IC.AItau,IC.time_constant, 'rows','complete');
subplot(1,3,3)
scatter(IC.AIrheobase,IC.rheobaseLP,5,'k')
xlabel('rheobase (AI)')
ylabel('rheobase')
[r_rheobase,p_rheobase] = corr(IC.AIrheobase,IC.rheobaseLP, 'rows','complete');
