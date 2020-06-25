%{
summaryICanalysis
%}

clear; close all; clc;                                                      % prepare Matlab workspace

% list of cells
mainFolder = '/home/dude/genpath/';                                                      % main folder for data (EDIT HERE)
cellList = dir([mainFolder, 'genpath/*.mat']);                                      % cell list
% free parameters (across sweep QC)
minNmaxThres = 10;                                                          % threshold for difference b/w min and max voltage
origStdThresMax = 4;                                                        % above this threshold we remove sweeps > 1.75 S.D.s
origStdThresMin = 0.8;                                                      % below this threshold cells are not QC'd sweep-wise any further

% summary file names
% savefilename = 'cell type details (2020 06 08)';                            % main data file name
% savefilenameInd = 'cell type details (2020 06 08) sp qc';                   % single cells file name (cell ID added below in loop)      

% load data
load('manual_entry_data.mat')                                               % load Michelle's NHP data
layerID = ID; dendrite_typeMJ = dendrite_type; access_resistance = ...      % adjust variable names for overlap
    acces_resistance;

clear ID dendrite_type acces_resistance                                     

load('cell_types_specimen_details.mat','donor__species','specimen__id',...
    'line_name','tag__dendrite_type','cell_reporter_status',...
    'structure__acronym','structure__layer');                               % load Allen Institute parameters

% initialize
qc_removed.List.Std = []; rmvdStdCount = 1;
qc_removed.List.MinMax = []; rmvdMMCount = 1;
qc_removed.List.NoRin = []; rmvdNoRinCount = 1;
qc_removed.List.NoRheo = []; rmvdNoRheoCount = 1;
qc_removed.List.InitRa = []; rmvdInitRaCount = 1;

qc.sp_mat = zeros(length(cellList),60);
qc.class_mat= zeros(length(cellList),75);
IC.input_current_s = zeros(length(cellList),75);


for n = 1:length(cellList)                                               % for each cells
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
        k  = find(specimen__id==str2num(cellID));

        % get Allen Institute parameters
        IC.ID(n,:) = ([num2str(specimen__id(k)),'__AI']);
        IC.specimen(n,1) = donor__species(k);
        IC.struct(n,1) = structure__acronym(k);
        IC.cortical_layer(n,1) = categorical(structure__layer(k));
        IC.transline(n,1) = line_name(k);
        IC.reporterStatus(n,1) = cell_reporter_status(k);
        IC.dendrite_type(n,1) = tag__dendrite_type(k);
        IC.access_resistance(n,1) = NaN;
        IC.temperature(n,1) = 0;
        
    % pull data for PCTD cell
    else                                                                    % if a JMT cell
        IC.ID(n,:) = cellID;
        IC.specimen(n,:) = categorical(cellstr('NHP'));
        IC.struct(n,:) = categorical(cellstr('PFC'));
        IC.transline(n,:) = categorical(cellstr('N/A'));
        IC.reporterStatus(n,:) = categorical(cellstr('N/A'));
        if sum(ismember(layerID,cellID))
            k = find(ismember(layerID,cellID)==1);
             if access_resistance(k,1) ~= 'N/A'                              
                 IC.access_resistance(n,1) = ...
                     str2num(char(access_resistance(k,1)));
             end
            IC.dendrite_type(n,1) = dendrite_typeMJ(k,1);
            IC.cortical_layer(n,1) = Layer(k,1);
             if temperature(k,1) ~= 'N/A'
                 IC.temperature(n,1) = str2num(char(temperature(k,1)));
             end
        else
            IC.access_resistance(n,1) = NaN;
            IC.dendrite_type(n,1) = categorical(cellstr('N/A'));
            IC.cortical_layer(n,1) = categorical(cellstr('N/A'));
            IC.temperature(n,1) = 0;
        end
    end

    % load our analysis
    load([mainFolder,cellList(n).name]);                         % load analysis file for cell
    
    % qc stuff
    qc_logic = zeros(1,8);                                                  % initialize QC matrix

    if a.LP.fullStruct == 1                                                 % if full data structure is available
        qc.ID{n,1} = cellID;                                                 % cell ID
        qc.V_vec(n,1:length(a.LP.rmp(1,:))) = round(a.LP.rmp(1,:),2);       % resting membrane potential
        qc.V_vecDelta(n,1:length(a.LP.rmp(1,:))) = ...
            round(a.LP.rmp(1,1) - a.LP.rmp(1,:),2);                         % diff RMP pre and post stimulus
        spqcmatn = zeros(length(a.LP.sweepAmps),10);                        % initialize count of QC removals matrix (each column is a criteria)
        binaryMatCount = 1;
        spqcvectag = nan(20,300);                                           % initialize QC tag storage
        input_current_spqc = zeros(20,1);                                   % initialize input current storage
        for k = 1:length(a.LP.sweepAmps)                                    % for each sweep
            
            % sweep-wise quality control parameters
            qc.restVpre(n,k) = round(a.LP.stats{k,1}.qc.restVPre,2);        % RMP pre stimulus
            qc.restVpost(n,k) = round(a.LP.stats{k,1}.qc.restVPost,2);      % RMP post stimulus
            qc.restVdiffpreNpost(n,k) = round( ...
                a.LP.stats{k,1}.qc.diffV_b_e,2);                            % diff RMP pre and post stimulus
            qc.rmse_pre_lt(n,k) = round(a.LP.stats{k,1}.qc.rmse_pre,2);     % long term RMS pre stimulus
            qc.rmse_post_lt(n,k) = round(a.LP.stats{k,1}.qc.rmse_post,2);   % lt RMS post stimulus
            qc.rmse_pre_st(n,k) = round(a.LP.stats{k,1}.qc.rmse_pre_st,2);  % st RMS pre stimulus
            qc.rmse_post_st(n,k) = round( ...
                a.LP.stats{k,1}.qc.rmse_post_st,2);                         % st RMS post stimulus
            qc_logic(1:6) = qc_logic(1:6)+a.LP.stats{k,1}.qc.logicVec;                % QC logic vector (each column is a criteria)
            
            % spike-wise QC processing            
           processSpQC                                                     % process spike-wise QC
            
            % assess the removal of this sweep
            if sum(a.LP.stats{k,1}.qc.logicVec) == 0 %isfield(a.LP,'stats') &&k<=length(a.LP.stats) && ...~isempty(a.LP.stats{k,1}) && ...
                qc.sweepID(n,k) = k;
                qc.sweepBinary(n,k) = 1;
                sweepBinaryOrig(1,k) = 1;
            else
                qc.sweepID(n,k) = 0;
                qc.sweepBinary(n,k) = 0;
                sweepBinaryOrig(1,k) = 0;
            end
            if isfield(a.LP.stats{k,1},'qcRemovals') && ...
                sum([sum(a.LP.stats{k,1}.qcRemovals.QCmatT2P), ...
                   sum(a.LP.stats{k,1}.qcRemovals.QCmatT2PRe), ...
                       sum(a.LP.stats{k,1}.qcRemovals.QCmatTrough), ...
                        sum(a.LP.stats{k,1}.qcRemovals.percentRheobaseHeight)]) > 0
              if 0.33*length(unique( ...
                      [ a.LP.stats{k, 1}.qcRemovals.minInterval,    ...
                        a.LP.stats{k, 1}.qcRemovals.dVdt0     ...
                        a.LP.stats{k, 1}.qcRemovals.mindVdt     ...   
                        a.LP.stats{k, 1}.qcRemovals.maxThreshold     ...
                        a.LP.stats{k, 1}.qcRemovals.minDiffThreshold2PeakN   ...
                        a.LP.stats{k, 1}.qcRemovals.minDiffThreshold2PeakB   ...
                        a.LP.stats{k, 1}.qcRemovals.diffthreshold2peakT   ...
                        a.LP.stats{k, 1}.qcRemovals.minIntervalRe   ...
                        a.LP.stats{k, 1}.qcRemovals.dVdt0Re   ...
                        a.LP.stats{k, 1}.qcRemovals.minTrough   ...
                        a.LP.stats{k, 1}.qcRemovals.percentRheobaseHeight   ...
                      ])) > length(a.LP.stats{k, 1}.spTimes) 
                  
                qc.sweeps_removed_SpQC(n,k) = k;
                qc.sweepBinary(n,k) = 0;
                qc.sweepID(n,k) = 0;
                qc.class_mat(n,k) = 10;
              end
            else
               qc.sweeps_removed_SpQC(n,k) = 0;
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
%             export_fig(['D:\my\genpath\',cellID,' spike QC (sweeps by critiera)'],'-pdf','-r100');
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
%             export_fig(['D:\my\genpath\',cellID,' spike QC (sweeps by spike binary)'],'-pdf','-r100');
%             close
%             clear spqcmatnbinaryid spqcmatnbinary k_len_spID
%         end

        qc.logic_mat(n,1:6) = qc_logic(1:6);
        if find(qc.sweeps_removed_SpQC)
         qc.logic_mat(n, 7) = length(find(qc.sweeps_removed_SpQC(n,:)));
        else
         qc.logic_mat(n, 7) = 0;
        end
        processBwSweepsQC                                                   % across sweep QC
    end
    if size(qc.class_mat,1)~=n
        qc.class_mat(n,:) = 0;
    end
    
    if size(qc.sweepBinary,1)==n && sum(qc.sweepBinary(n,:))>0
        
        % subthreshold summary parameters
        IC.acquireRes(n,1) = double(a.LP.acquireRes);
        IC.resistance(n,1) = round(double(a.LP.subSummary.resistance),2);            % Allen Method
        [IC.resistance_hd(n,1), IC.rect_index(n,1)]  = resistance_hd(a.LP, cellID);       % resistance (HD?)
        IC.resistance_ss(n,1) = resistance_ss(a.LP, cellID);                        % resistance based on steady state
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
                    if IC.slowTrough(n,1) < 2                               %  check if the cell has a true fAHP, its peak should be more depolarized than any other afterhyperpolarizations                         
                        IC.fAHPamp(n,1) = IC.thresholdLP(n,1) - ...
                            IC.trough(n,1);                                 %  amplitude of fAHP from threshold
                    else 
                        IC.fAHPamp(n,1) = NaN;
                    end
                    if size(a.LP.stats{idx(k),1}.waves,1) > 1
                        wf(n,:) = round(mean(a.LP.stats{idx(k),1}.waves),2);
                    else
                        wf(n,:) = round(a.LP.stats{idx(k),1}.waves,2);
                    end
                    break
                elseIC
                    IC.rheobaseLP(n,1) = NaN;
                end
            end
            
            % global spike parameters and Hero sweep selection
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
                    IC.rate_1s_hero(n,1) = NaN;
                    IC.burst_hero(n,1) = NaN;
                    IC.delay_hero(n,1) =   NaN;
                    IC.latency_hero(n,1) = NaN;
                    IC.cv_ISI(n,1) =  NaN;
                    IC.adaptation2(n,1) = NaN;
                    IC.hero_amp(n,1) = NaN;
                elseif length(k) > 1
                    IC.rate_1s_hero(n,1) = mean(IC.rate_1s(n,find(IC.input_current_s(n,:)==a.LP.sweepAmps(k))));
                    IC.burst_hero(n,1) = mean(train_burst(n,k(1:length(k))));
                    IC.delay_hero(n,1) =   mean(train_delay(n,k(1:length(k))));
                    IC.latency_hero(n,1) = mean(train_latency(n,k(1:length(k))));
                    IC.cv_ISI(n,1) =   mean(train_cv_ISI(n,k(1:length(k))));
                    IC.adaptation1(n,1) = mean(train_adaptation1(n,k(1:length(k))));
                    IC.adaptation2(n,1) = mean(train_adaptation2(n,k(1:length(k))));
                    IC.hero_amp(n,1) = unique(a.LP.sweepAmps(k));       
                else
                    IC.rate_1s_hero(n,1) = IC.rate_1s(n,find(IC.input_current_s(n,:)==a.LP.sweepAmps(k)));
                    IC.burst_hero(n,1) =  a.LP.stats{k, 1}.burst;
                    IC.delay_hero(n,1) =   a.LP.stats{k, 1}.delay;
                    IC.latency_hero(n,1) = a.LP.stats{k, 1}.latency;
                    IC.cv_ISI(n,1) =   a.LP.stats{k, 1}.cvISI ;  
                    IC.adaptation1(n,1) = a.LP.stats{k, 1}.adaptIndex; 
                    IC.adaptation2(n,1) = a.LP.stats{k, 1}.adaptIndex2; 
                    IC.hero_amp(n,1) = a.LP.sweepAmps(k);
                end    
                clear B I idx amp temp k
            end
        end
    end
end
clear cell_reporter_status donor__species line_name specimen__id structure__acronym ...
    structure__layer tag__dendrite_type

% Cleaning up variables 
fieldnames_var = fieldnames(IC);   

for n = 1:length(cellList) 
    for var = 1:length(fieldnames_var)  
            if isnumeric(IC.(fieldnames_var{var})(n,1)) && IC.(fieldnames_var{var})(n,1) == 0
               IC.(fieldnames_var{var})(n,1) = NaN;
            end
    end
end
% Getting the variable names to overwrite in them

for n = 1:length(cellList)      
 if isnan(IC.resistance_ss(n,1)) && isnan(IC.resistance_ss(n,1))
    for var = 11:length(fieldnames_var)                                % Leave the first 12 variables untouched, since they are still usefull
         IC.(fieldnames_var{var})(n,1) = NaN;
    end
   qc_removed.List.NoRin{rmvdNoRinCount,1} = IC.ID(n,:); 
   rmvdNoRinCount = rmvdNoRinCount + 1;    
 end
 if  isnan(IC.rheobaseLP(n,1))
    for var = 11:length(fieldnames_var)                                % Leave the first 12 variables untouched, since they are still usefull
         IC.(fieldnames_var{var})(n,1) = NaN;
    end
     IC.(fieldnames_var{var})(n,1) = NaN;
     qc_removed.List.NoRheo{rmvdNoRheoCount,1} = IC.ID(n,:); 
     rmvdNoRheoCount = rmvdNoRheoCount + 1;
 end
 if IC.access_resistance(n,1) > 20    
    for var = 11:length(fieldnames_var)                                % Leave the first 12 variables untouched, since they are still usefull
            IC.(fieldnames_var{var})(n,1) = NaN;     
    end
   qc.sweepBinary(n,:) = 0;  
   qc_removed.List.InitRa{rmvdInitRaCount,1} = IC.ID(n,:); 
   rmvdInitRaCount = rmvdInitRaCount + 1;    
 end   
end

variable_range = [10:length(fieldnames_var)];                              %selection of variables fron Eric IC for new IC
memory = [];
IC.ID = string(IC.ID);


for  n = 1:length(IC.ID)                                                   % these two loops determine which cell ID should be removed
 check = 0 ;                                                               % resets the check variable which counts the numbers of NaNs in variables
 for var = variable_range(1,:)
    check = check + isnan(IC.(fieldnames_var{var})(n,1));                  % goes through the variables and adds a one to check if there is a NaN
 end
 if check >= 60 %length(variable_range(1,:))                                   % if there are more NaNs then non-manual-data-entry variables
    memory = [memory , n];                                                 % remember the cell
 end
end

qc_removed.IDs =  IC.ID(memory,:);
qc_removed.sp_mat = qc.sp_mat(memory,:);
qc_removed.logic_mat = qc.logic_mat(memory,:); 
qc_removed.sweepBinary = qc.sweepBinary(memory,:);  
qc_removed.rmse_pre_lt = qc.rmse_pre_lt(memory,:); 
qc_removed.rmse_post_lt = qc.rmse_post_lt(memory,:); 
qc_removed.restVdiffpreNpost = qc.restVdiffpreNpost(memory,:);  
qc_removed.sweeps_removed_SpQC = qc.sweeps_removed_SpQC(memory,:);
qc_removed.List.RMS_exHigher9 = qc_removed.IDs(max(qc_removed.logic_mat(:,3:4),[], 2) > 9);
qc_removed.List.Depol_exHigher8 = qc_removed.IDs(qc_removed.logic_mat(:,6) > 8);
qc_removed.List.justNoRheo = setdiff(qc_removed.List.NoRheo, qc_removed.List.NoRin);
qc_removed.List.justNoRin = setdiff(qc_removed.List.NoRin, qc_removed.List.NoRheo);

qc = structfun(@(x) (removerows(x, 'ind', memory)), qc, 'UniformOutput', false);
IC = structfun(@(x) (removerows(x, 'ind', memory)), IC, 'UniformOutput', false); % delete out all cells indicated by the variable memory 
 
clearvars -except qc qc_removed IC  

return

