%{
summaryICanalysis
%}

clear; close all; clc;                                                      % prepare Matlab workspace

% list of cells
mainFolder = 'D:\patchseq_matlab\';                                                      % main folder for data (EDIT HERE)
cellList = dir([mainFolder,'\*.mat']);                                      % cell list
    
% free parameters (across sweep QC)
minNmaxThres = 10;                                                          % threshold for difference b/w min and max voltage
origStdThresMax = 3;                                                        % above this threshold we remove sweeps > 1.75 S.D.s
origStdThresMin = 0.7;                                                      % below this threshold cells are not QC'd sweep-wise any further

% load data
load('mouse metadata.mat','cell_specimen_id','ephys_session_id', ...
    'dendrite_type','depth_from_pia_um','structure','cell_specimen_name');  % load Allen Institute parameters

% initialize
IC.input_current_s = zeros(length(cellList),100);
IC.firing_rate_s = zeros(length(cellList),100);
removedListStd = []; rmvdStdCount = 1;
removedListMinMax = []; rmvdMMCount = 1;
spqcmat = zeros(length(cellList),100);

for n = 1:length(cellList)                                                  % for each cells
    clc; disp(n)                                                            % display n value
    sweepIDcount = 1;
    
    % load our analysis
    load([mainFolder,'\',cellList(n).name]);                            % load analysis file for cell
    
    % match species ID # to get Allen parameters
    if n < 6
        ephysID = cellList(n).name(length(cellList(n).name)-21:length(cellList(n).name)-12);
    else
        ephysID = cellList(n).name(length(cellList(n).name)-20:length(cellList(n).name)-12);
    end
    k  = find(ephys_session_id==str2num(ephysID));

    % get Allen Institute parameters
    IC.cellID{n,:} = num2str(cell_specimen_id(k));
    IC.ephysID{n,:} = ephysID;
    IC.struct(n,1) = structure(k);
    IC.cort_depth(n,1) = depth_from_pia_um(k);
    temp = convertStringsToChars(cell_specimen_name(k));
    IC.transline(n,1) = convertCharsToStrings(temp(1:3));
    IC.dendrite_type(n,1) = dendrite_type(k);

    if a.LP.fullStruct==1 && n ~= 2013 %more sweepAmps than restVPre

        % qc stuff
        qc_logic = zeros(1,9);                                              % initialize QC matrix

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
            processSpQC                                                     % process spike-wise QC
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
%         if sum(sum(spqcmatn))>0                                             % figure of spike-wise QC (sweeps X criteria)
%             figure('Position',[50 50 300 300]); set(gcf,'color','w');
%             imagesc(spqcmatn)
%             box off
%             colorbar
%             colormap('gray')
%             xticks(1:10)
%             xticklabels({'interval','null dV/dt','dV/dt<5mV/ms', ...
%                 'threshold>-20mV','t2pN(B)<35(45)mV','t2pT>1.5ms','interval Re', ...
%                 'null dV/dt Re','trough>-30mV','<30% Rheobase height'})
%             xtickangle(45)
%             ylabel('current input (pA)')
%             yticks(1:size(spqcmatn,1))
%             yticklabels({a.LP.sweepAmps})
%     %             export_fig(['D:\my\genpath\',cellID,' spike QC (sweeps by critiera)'],'-pdf','-r100');
%             close
%         end 
%         if exist('spqcmatnbinaryid','var')                                  % figure of spike-wise QC (sweeps X spikes)
%             figure('Position',[50 50 300 250]); set(gcf,'color','w');
%             imagesc(spqcmatnbinaryid)
%             box off
%             colormap('gray')
%             xlabel('spike # (white==passes QC)')
%             ylabel('current input (pA)')
%             yticks(1:length(k_len_spID))
%             yticklabels({k_len_spID})
%     %             export_fig(['D:\my\genpath\',cellID,' spike QC (sweeps by spike binary)'],'-pdf','-r100');
%             close
            clear spqcmatnbinaryid spqcmatnbinary k_len_spID
%         end
        % add sp qc
        qc_logic_mat(n,1:9) = qc_logic;
        processBwSweepsQC                                                   % across sweep QC

        if size(sweepBinary,1)==n && sum(sweepBinary(n,:))>0

            % subthreshold summary parameters
            IC.acquireRes(n,1) = mean(double(a.LP.acquireRes));
            IC.acquireRes(n,2) = std(double(a.LP.acquireRes));
            IC.resistance(n,1) = round(double(a.LP.subSummary.resistance),2);   % Allen Method
            IC.resistance_hd(n,1) = resistance_hd(a.LP);                        % resistance (HD?)
            IC.resistance_ss(n,1) = resistance_ss(a.LP);                        % resistance based on steady state
            IC.Vrest(n,1) = restingMP(a.LP);                                    % resting membrane potential
            IC.time_constant(n,1) = round(double(a.LP.subSummary.tauMin),2);
            IC.time_constant2(n,1) = round(double(a.LP.subSummary.tauSS),2);

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
    %                 IC.rate_750ms(n,k) = a.LP.stats{idx(k),1}.meanFR750;
    %                 IC.rate_500ms(n,k) = a.LP.stats{idx(k),1}.meanFR500;
    %                 IC.rate_250ms(n,k) = a.LP.stats{idx(k),1}.meanFR250;
    %                 IC.rate_100ms(n,k) = a.LP.stats{idx(k),1}.meanFR100;
    %                 IC.rate_50ms(n,k) = a.LP.stats{idx(k),1}.meanFR50;
                    IC.delay(n,k) = round(double(a.LP.stats{idx(k),1}.delay(1,1)),2);
                    IC.burst(n,k) = round(double(a.LP.stats{idx(k),1}.burst(1,1)),2);
                    IC.latency(n,k) = round(double(a.LP.stats{idx(k),1}.latency(1,1)),2);
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
    %                 IC.rate_750msRheobase(n,1) = a.LP.stats{idx(k),1}.meanFR750;
    %                 IC.rate_500msRheobase(n,1) = a.LP.stats{idx(k),1}.meanFR500;
    %                 IC.rate_250msRheobase(n,1) = a.LP.stats{idx(k),1}.meanFR250;
    %                 IC.rate_100msRheobase(n,1) = a.LP.stats{idx(k),1}.meanFR100;
    %                 IC.rate_50msRheobase(n,1) = a.LP.stats{idx(k),1}.meanFR50;
                    IC.delayRheobase(n,1) = round(double(a.LP.stats{idx(k),1}.delay(1,1)),2);
                    IC.burstRheobase(n,1) = round(double(a.LP.stats{idx(k),1}.burst),2);
                    IC.latencyRheobase(n,1) = round(double(a.LP.stats{idx(k),1}.latency(1,1)),2);
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
%                     figure('Position',[50 50 300 250]); set(gcf,'color','w');
%                     plot(a.LP.V{1,idx(k)},'k')
%                     hold on
%                     scatter(IC.peakTimeLP(n,1),IC.peakLP(n,1),30,'g')
%                     scatter(IC.peakTimeLP(n,1)+(IC.fastTrough(n,1)/IC.acquireRes(n,1)),a.LP.stats{idx(k),1}.fastTrough(1),30,'r')
%                     scatter(IC.peakTimeLP(n,1)+(IC.slowTrough(n,1)/IC.acquireRes(n,1)),a.LP.stats{idx(k),1}.slowTrough(1),50,'k')
%                     xlabel('time')
%                     ylabel('voltage (mV)')
%                     axis tight
%                     xlim([a.LP.stimOn(1,idx(k))-(1/IC.acquireRes(n,1)) ...
%                         IC.peakTimeLP(n,1)+((IC.slowTrough(n,1)+1)/IC.acquireRes(n,1))])
%                     box off
%                     annotation('textbox',[.3 .7 0 0],'String',IC.fastTrough(n,1),'FitBoxToText','on','LineStyle','none','color','r');
%                     annotation('textbox',[.3 .9 0 0],'String',IC.slowTrough(n,1),'FitBoxToText','on','LineStyle','none');
%                     annotation('textbox',[.52 .7 0 0],'String',IC.fastTroughV(n,1),'FitBoxToText','on','LineStyle','none','color','r');
%                     annotation('textbox',[.52 .9 0 0],'String',IC.slowTroughV(n,1),'FitBoxToText','on','LineStyle','none');
%                     export_fig(['D:\figs\',IC.cellID{n,:},' ',int2str(idx(k)),' trough'],'-pdf','-r300');
%                     close

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
                IC.mdn_insta_freq(n,1) = median_isi(a.LP);                  % obtain the median ISI of all suprathreshold sweeps

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
    %                 IC.rate_750msHero(n,1) = mean(a.LP.stats{k,1}.meanFR750);
    %                 IC.rate_500msHero(n,1) = mean(a.LP.stats{k,1}.meanFR500);
    %                 IC.rate_250msHero(n,1) = mean(a.LP.stats{k,1}.meanFR250);
    %                 IC.rate_100msHero(n,1) = mean(a.LP.stats{k,1}.meanFR100);
    %                 IC.rate_50msHero(n,1) = mean(a.LP.stats{k,1}.meanFR50);
    %                 IC.rateIntWindows(n,:) = a.LP.stats{k,1}.frIntWindows;
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
    %                 IC.rate_750msHero(n,1) = a.LP.stats{k,1}.meanFR750;
    %                 IC.rate_500msHero(n,1) = a.LP.stats{k,1}.meanFR500;
    %                 IC.rate_250msHero(n,1) = a.LP.stats{k,1}.meanFR250;
    %                 IC.rate_100msHero(n,1) = a.LP.stats{k,1}.meanFR100;
    %                 IC.rate_50msHero(n,1) = a.LP.stats{k,1}.meanFR50;
    %                 IC.rateIntWindows(n,:) = a.LP.stats{k,1}.frIntWindows;
                    IC.burst_hero(n,1) =  a.LP.stats{k,1}.burst;
                    IC.delay_hero(n,1) =   a.LP.stats{k,1}.delay(1,1);
                    IC.latency_hero(n,1) = a.LP.stats{k,1}.latency(1,1);
                    IC.cv_ISI(n,1) =   a.LP.stats{k,1}.cvISI ;  
                    IC.adaptation1(n,1) = a.LP.stats{k,1}.adaptIndex; 
                    IC.adaptation2(n,1) = a.LP.stats{k,1}.adaptIndex2; 
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
