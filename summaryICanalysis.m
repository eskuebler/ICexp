%{
summaryICanalysis
%}

clear; close all; clc;                                                      % prepare Matlab workspace

% list of cells
mainFolder = 'D:\my\';                                                      % main folder for data (EDIT HERE)
cellList = dir([mainFolder,'genpath\*.mat']);                                 % cell list
    
% free parameters (across sweep QC)
minNmaxThres = 10;                                                          % threshold for difference b/w min and max voltage
origStdThresMax = 3;                                                        % above this threshold we remove sweeps > 1.75 S.D.s
origStdThresMin = 0.7;                                                      % below this threshold cells are not QC'd sweep-wise any further      

% load data
load('manual_entry_data.mat')                                               % load Michelle's NHP data
layerID = ID; dendrite_typeMJ = dendrite_type; clear ID dendrite_type       % adjust variable names for overlap
load('pyramidal cell IDs.mat')                                              % load Michelle's pyramidal cells
pyrID = ID; clear ID                                                        % adjust variable names for overlap
load('cell_types_specimen_details.mat','donor__species','specimen__id',...
    'line_name','tag__dendrite_type',...
    'ef__threshold_i_long_square','cell_reporter_status',...
    'structure__acronym','structure__layer');                               % load Allen Institute parameters

% initialize
IC.input_current_s = zeros(length(cellList),75);
IC.firing_rate_s = zeros(length(cellList),75);
removedListStd = []; rmvdStdCount = 1;
removedListMinMax = []; rmvdMMCount = 1;
spqcmat = zeros(length(cellList),60);

for n = 1:length(cellList)                                                  % for each cells
    clc; disp(n)                                                            % display n value
    sweepIDcount = 1;
    
    %% Determining the from which institution the cell was obtained and assigning their matadata and control parameters
    
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
        IC.ID(n,:) = {([num2str(specimen__id(k)),'__AI'])};
        IC.specimen(n,1) = donor__species(k);
        IC.struct(n,1) = structure__acronym(k);
        IC.cortical_layer(n,1) = categorical(structure__layer(k));
        IC.transline(n,1) = line_name(k);% Michelle Changes
        IC.reporterStatus(n,1) = cell_reporter_status(k);
        IC.dendrite_type(n,1) = tag__dendrite_type(k);
        AIrheobase(n,1) = round(ef__threshold_i_long_square(k));                                                                                     % Michelle change; deleted IC.PyrID
        IC.access_resistance(n,1) = NaN;
        IC.temperature(n,1) = 0;
        
    % pull data for PCTD cell
    else                                                                    % if a JMT cell
        IC.ID(n,:) = {cellID} ;
        IC.specimen(n,:) = categorical(cellstr('NHP'));
        IC.struct(n,:) = categorical(cellstr('PFC'));
        IC.transline(n,:) = categorical(cellstr('N/A'));
        IC.reporterStatus(n,:) = categorical(cellstr('N/A'));
        
        if sum(ismember(layerID,cellID))
            k = find(ismember(layerID,cellID)==1);
            
            %Not implemented, because values have to be revised
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
   
%% QC filters  

    % load our analysis
    load([mainFolder,'genpath\',cellList(n).name]);                         % load analysis file for cell
    
    % qc stuff
    qc_logic = zeros(1,6);                                                  % initialize QC matrix

    if a.LP.fullStruct == 1                                                 % if full data structure is available
        qcID{n,1} = cellID;                                                 % cell ID
        qc_V_vec(n,1:length(a.LP.rmp(1,:))) = round(a.LP.rmp(1,:),2);       % resting membrane potential
       % spqcmatn = zeros(length(a.LP.sweepAmps),6);                         % initialize count of QC removals matrix (each column is a criteria)
       % spqcmatnbinary = nan(20,300); spqcmatnbinaryid = nan(20,300);       % initialize spike-wise QC matrix, assuming that there are no more than 300 spikes per sweep
       % spqcvectag = nan(20,300);                                           % initialize QC tag storage
       %  input_current_spqc = zeros(20,1);                                   % initialize input current storage
        for k = 1:length(a.LP.sweepAmps)                                    % for each sweep
            
            % spike-wise QC processing            
            %processSpQC                                                     % process spike-wise QC
            % remove sweeps that exceed good/bad spike ratio  0.3
            % add to QC vec                  
           
            % assess the removal of this sweep
            if sum(a.LP.stats{k,1}.qc.logicVec) == 0 %isfield(a.LP,'stats') &&k<=length(a.LP.stats) && ...~isempty(a.LP.stats{k,1}) && ...
                sweepID(n,k) = k;
                sweepBinary(n,k) = 1;
                sweepBinaryOrig(1,k) = 1;
            else
                sweepID(n,k) = 0;
                sweepBinary(n,k) = 0;
                sweepBinaryOrig(1,k) = 0;
            end
        end
        
        % add sp qc
        qc_logic_mat(n,1:6) = qc_logic;
        processBwSweepsQC                                                   % across sweep QC
    end
    if size(qc_class_mat,1)~=n
         qc_class_mat(n,:) = 0;
    end
    
%     if exist('input_current_spqc','var')                                    % save csv files for spike-wise QC
%         T = table(input_current_spqc,spqcmatnbinary);
%         writetable(T,[mainFolder,'genpath\',savefilenameInd,' ',cellID,...
%             '.xlsx'],'Sheet','Sheet1','WriteRowNames',true)
%         T = table(input_current_spqc,spqcmatnbinaryid);
%         writetable(T,[mainFolder,'genpath\',savefilenameInd,' ',cellID,...
%             '.xlsx'],'Sheet','Sheet2','WriteRowNames',true)
%         T = table(input_current_spqc,spqcvectag);
%         writetable(T,[mainFolder,'genpath\',savefilenameInd,' ',cellID, ...
%             '.xlsx'],'Sheet','Sheet3','WriteRowNames',true)
%        % T = table(input_current_spqc,spqcmatn);
%        % writetable(T,[mainFolder,'genpath\',savefilenameInd,' ',cellID, ...
%           %  '.xlsx'],'Sheet','Sheet4','WriteRowNames',true)
%         clear input_current_spqc spqcmatnbinary spqcmatnbinaryid
%     end
    
%    if size(sweepBinary,1)==n && sum(sweepBinary(n,:))>0
                
%         if n == 294
%            'YO' 
%         end

%% subthreshold summary parameters
IC.resistance_hd(n,1) = Michelle_calc_resistance_hd(a.LP); % Michelle Changes new line calling new function
IC.resistance_ss(n,1) = Michelle_calc_resistance_ss(a.LP); % Michelle Changes new line calling new function
IC.Vrest(n,1) = Michelle_calc_rmp(a.LP);                    % Michelle Changes new line calling new function

      if a.LP.fullStruct == 1          
        acquireRes(n,1) = double(a.LP.acquireRes);
   
        IC.time_constant(n,1) = round(double(a.LP.subSummary.tauMin),2);% Michelle Changes
        %time_constant2(n,1) = round(double(a.LP.subSummary.tauSS),2);
        
            k = find(ismember(sweepID(n,:),...
                find(round(double(a.LP.sweepAmps)) == -90))==1);            % find -90 pA input
            if length(k)>1
                k = k(1);
            end
            if ~isempty(k)
                getSubthresholdStats                                        % get subthreshold stats
            else                                                            % if no -90 pA sweep
                k = find(ismember(sweepID(n,:),...
                    find(round(double(a.LP.sweepAmps)) == -70))==1);       % find -70 pA sweep % Michelle Changes
                if length(k)>1
                    k = k(1);
                end
                if ~isempty(k)
                    getSubthresholdStats                                    % get subthreshold stats
                else                                                        % if no -70 pA sweep
                    k = find(ismember(sweepID(n,:),...
                        find(round(double(a.LP.sweepAmps)) == -110))==1);    % find -110 pA sweep 
                    if length(k)>1
                        k = k(1);
                    end
                    if ~isempty(k)
                        getSubthresholdStats                                % get subthreshold stats
                    else                                                    % if no -110pA sweeps
                        IC.subamp(n,1) = NaN;                                  % add NaNs (blank spaces in csv format)
                        IC.submin(n,1) = NaN;
                       % rebound_slope(n,1) = NaN;
                       % rebound_depolarization(n,1) = NaN;
                       % nb_rebound_sp(n,1) = 0;
                        IC.sag(n,1) = NaN;
                        IC.steadystate(n,1) = NaN;
                        IC.sag_ratio(n,1) = NaN;
                    end
                end
            end
      end
             
            %% find rheobase and parameters of first spike
            
      
      if a.LP.fullStruct == 1      
            [B,I] = sort(round(double(a.LP.sweepAmps(1:length(a.LP.stats)))));
            input_current_all(n,1:length(B)) = B;
            int_vec = find(B>0);
            temp = int_vec(find(ismember(int_vec,sweepID(n,:))==1));
            idx = I(temp);
            amp = B(temp);
            IC.input_current_s(n,1:length(amp)) = round(double(a.LP.sweepAmps(idx)));

            for k = 1:length(idx)
                if isfield(a.LP.stats{idx(k),1},'spTimes') && ...
                        sum(~isnan(a.LP.stats{idx(k),1}.spTimes))>0
                    
                    % spike train parameters
                    IC.firing_rate_s(n,k) = a.LP.stats{idx(k),1}.meanFR1000;
                    train_delay(n,k) = round(double(a.LP.stats{idx(k),1}.delay),2);
                    train_burst(n,k) = round(double(a.LP.stats{idx(k),1}.burst),2);
                    train_latency(n,k) = round(double(a.LP.stats{idx(k),1}.latency),2);
                    train_cv_ISI(n,k) = round(double(a.LP.stats{idx(k),1}.cvISI),2);
                    train_adaptation1(n,k) = round(double(a.LP.stats{idx(k),1}.adaptIndex),2);
                    train_adaptation2(n,k) = round(double(a.LP.stats{idx(k),1}.adaptIndex2),2);
                    train_peak_adaptation1(n,k) = round(double(a.LP.stats{idx(k),1}.peakAdapt),2);
                    train_peak_adaptation2(n,k) = round(double(a.LP.stats{idx(k),1}.peakAdapt2),2);
                    peak_adaptation1(n,1) = round(double(a.LP.stats{idx(k),1}.peakAdapt),2);
                    peak_adaptation2(n,1) = round(double(a.LP.stats{idx(k),1}.peakAdapt2),2);
                    
                    
                end
            end
            
             for k = 1:length(idx)
                if isfield(a.LP.stats{idx(k),1},'spTimes') && ...
                        sum(~isnan(a.LP.stats{idx(k),1}.spTimes))>0
                    % spike analysis
                    IC.rheobaseLP(n,1) = amp(k);
                    IC.delay_rheo(n,1) = round(double(a.LP.stats{idx(k),1}.delay),2);
                    IC.burst_rheo(n,1) = round(double(a.LP.stats{idx(k),1}.burst),2);
                    IC.latency_rheo(n,1) = round(double(a.LP.stats{idx(k),1}.latency),2);
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
                    if IC.slowTrough(n,1) < 2                                %  check if the cell has a true fAHP, its peak should be more depolarized than any other afterhyperpolarizations                         
                     IC.fAHPamp(n,1) = IC.thresholdLP(n,1) - IC.trough(n,1); %  amplitude of fAHP from threshold
                    else 
                     IC.fAHPamp(n,1) = NaN;
                    end
%                    if size(a.LP.stats{idx(k),1}.waves,1) > 1
%                        wf(n,:) = round(mean(a.LP.stats{idx(k),1}.waves),2);
%                    else
%                        wf(n,:) = round(a.LP.stats{idx(k),1}.waves,2);
%                    end
%                    break
                else
                    IC.rheobaseLP(n,1) = NaN;
                end
            end
      end
%% Global spiketrain parameters and picking "Hero sweep" for more spike train parameters per cell  
              
     if a.LP.fullStruct == 1 
        IC.maxFiringRate(n,1) = max(IC.firing_rate_s(n,:));                % Obtain maximum firing rate
        IC.mdn_insta_freq(n,1) = Michelle_median_isi(a.LP);                % Obtain the median ISI of all suprathreshold sweeps
        k = [];                                                            % resetting k
        flag = 0;                                                          % variable to fire the if condition in while loop only one time
        
        if  length(IC.rheobaseLP) == n 
        [~, k] = min(abs(double(B)-IC.rheobaseLP(n,1)*1.5));  
           if k > 0  
             while sum(ismember(sweepID(n,:), k)) == 0 || ...
                     isfield(a.LP.stats{k,1},'burst') == 0 || ...
                     IC.firing_rate_s(n,k)== 0 || ...
                     a.LP.sweepAmps(k) <= IC.rheobaseLP(n,1) || ...
                     a.LP.sweepAmps(k) > 2*IC.rheobaseLP(n,1)
                  k = k -1;  
                  if k == 0 && flag == 0
                   [~, k] = min(abs(double(B)-IC.rheobaseLP(n,1)*8)); 
                   flag = 1;
                  end
                  if k== 0 && flag == 1
                      break
                  end
             end  
           end
           if k == 0 
                  IC.firing_rate_s_hero(n,1) = NaN;
                  IC.burst_hero(n,1) = NaN;
                  IC.delay_hero(n,1) =   NaN;
                  IC.latency_hero(n,1) = NaN;
                  IC.cv_ISI(n,1) =  NaN;
                  IC.adaptation2(n,1) = NaN;
                  IC.hero_amp(n,1) = NaN;
           elseif length(k) > 1
           IC.firing_rate_s_hero(n,1) = mean(IC.firing_rate_s(n,k(1:length(k))));
           IC.burst_hero(n,1) = mean(train_burst(n,k(1:length(k))));
           IC.delay_hero(n,1) =   mean(train_delay(n,k(1:length(k))));
           IC.latency_hero(n,1) = mean(train_latency(n,k(1:length(k))));
           IC.cv_ISI(n,1) =   mean(train_cv_ISI(n,k(1:length(k))));
           IC.adaptation1(n,1) = mean(train_adaptation1(n,k(1:length(k))));
           IC.adaptation2(n,1) = mean(train_adaptation2(n,k(1:length(k))));
           IC.hero_amp(n,1) = unique(a.LP.sweepAmps(k));       
           else
           IC.firing_rate_s_hero(n,1) = IC.firing_rate_s(n,k(1:length(k)));
           IC.burst_hero(n,1) =  a.LP.stats{k, 1}.burst;
           IC.delay_hero(n,1) =   a.LP.stats{k, 1}.delay;
           IC.latency_hero(n,1) = a.LP.stats{k, 1}.latency;
           IC.cv_ISI(n,1) =   a.LP.stats{k, 1}.cvISI ;  
           IC.adaptation1(n,1) = a.LP.stats{k, 1}.adaptIndex; 
           IC.adaptation2(n,1) = a.LP.stats{k, 1}.adaptIndex2; 
           IC.hero_amp(n,1) = a.LP.sweepAmps(k);
           
           end         
        end
     end 
     
end

%% Cleaning up variables 

IC.resistance_ss(IC.resistance_ss==0)= NaN;                                % 0 is a non-sensical value in most variables and should be NaN
IC.resistance_hd(IC.resistance_hd==0)= NaN;
IC.time_constant(IC.time_constant==0)= NaN;
IC.Vrest_sag_sweep(IC.Vrest_sag_sweep==0)= NaN;                           
IC.sag_ratio(IC.sag_ratio==0)= NaN;                                        
IC.mdn_insta_freq(IC.mdn_insta_freq==0)= NaN;  
IC.rheobaseLP(IC.rheobaseLP==0)= NaN;  
IC.hero_amp(IC.hero_amp==0)=NaN;
IC.firing_rate_s(IC.firing_rate_s==0)= NaN;

fieldnames_var = fieldnames(IC);                                           % Getting the variable names to overwrite in them

for  n = 1:length(cellList)  
  if isnan(IC.resistance_hd(n,1)) || isnan(IC.rheobaseLP(n,1))             % If crucial features cannot be determined, all parameters are set to NaN  
    for var = 10:length(fieldnames_var)-2                                   % Leave the first 7 and last 2 untouched, since they are still usefull
      IC.(fieldnames_var{var})(n,1) = NaN;
    end
  end
 end



