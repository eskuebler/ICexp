%{

mainCatalogJMT.m - intracellular experiment catalog

*** this script highlights a huge issue with incorporating data across labs
recorders have made many different names for the protocols, which are 
used to ID what experiment was done, going forward we need a set list and 
set names that do not get adjusted, or at the very least get noted as 
adjusted ***

%}

clear; close all; clc;                                                      % prepare matlab 

files = dir('H:\NHP cell type database\**\*.abf');                          % get file info

% initialize storage and protocol counts
tstamp = zeros(length(files),1);                                            % time-stamps for start of recording
protocolName = cell(length(files),1);                                       % protocol name
LPcount = 1; SPcount = 1; GFcount = 1; SSPcount = 1; RampUpcount = 1;
NAcount = 1; RampNewcount = 1; Noisecount = 1; NotAvailcount = 1; 
CCcount = 1; VCcount = 1;

% for each recording file
for k = 1:length(files)                                                     % for each file
    display(files(k,1).name)                                                % display
    file = [files(k,1).folder,'\',files(k,1).name];                         % set folder/file name
    [~,~,h] = abfload(file,'sweeps','a','channels','a');                    % load ABF data
    tstamp(k,1) = h.uFileStartDate;                                         % set time of experiment start

    % statement to determine which protocol was used (many different names
    % here for some of the protocols (i.e., 1000 ms, 3 ms)
    if sum(h.protocolName(end-22:end-4) == 'Monkey_1000 ms step') == 19 ...
            || length(h.protocolName)>=28 ...
            && sum(h.protocolName(end-26:end-4) == 'Monkey_1000 ms steplong')==23 ...
            || length(h.protocolName)>=28 ...
            && sum(h.protocolName(end-27:end-4) == 'Monkey_1000 ms stepExtra')==24
        
        protocolName{k,1} = '1000 ms step';
        LPcount = LPcount + 1;                                              % count of long pulse abf
    
    elseif sum(h.protocolName(end-19:end-4) == 'Monkey_3 ms step')==16 ...
            || sum(h.protocolName(end-20:end-4) == 'Monkey_3 ms step2')==17 ...
            || length(h.protocolName)>=28 ...
            && sum(h.protocolName(end-28:end-4) == 'Monkey_3 ms step_addsweep')==25
        
        protocolName{k,1} = '3 ms step';
        SPcount = SPcount + 1;                                              % count of short pulse abf
    
    elseif sum(h.protocolName(end-10:end-4) == 'Gapfree') == 7 ...
            || sum(h.protocolName(end-11:end-4) == 'Gap free') == 8 ...
            || sum(h.protocolName(end-18:end-4) == 'Monkey Gap free') == 15
    
        protocolName{k,1} = 'Gapfree';
        GFcount = GFcount + 1;                                              % count of gap free
    
    elseif length(h.protocolName)>=30 ...
            && sum(h.protocolName(end-30:end-4) == 'Monkey_0.5ms square pulse 1')==27
        
        protocolName{k,1} = '0.5ms square pulse';
        SSPcount = SSPcount + 1;                                            % count of 0.5 ms pulse
        
    elseif length(h.protocolName)>=28 ...
            && sum(h.protocolName(end-28:end-4) == 'Monkey_25 pA per sec ramp')==25
    
        protocolName{k,1} = '25 pA per sec ramp';
        RampUpcount = RampUpcount + 1;                                      % count ramp pulse (AI)
        
    elseif sum(h.protocolName(end-21:end-4) == 'Monkey_CC_RAMP_NEW')==18
    
        protocolName{k,1} = 'ramp new';
        RampNewcount = RampNewcount + 1;                                    % count ramp new (EPSP) 
        
    elseif sum(h.protocolName(end-19:end-4) == 'Monkey_1nA_pulse')==16
    
        protocolName{k,1} = '1nA_pulse';                                    % count nanoamp pulse
        NAcount = NAcount + 1; 
        
    elseif sum(h.protocolName(end-21:end-4) == 'CC_Capacitance1sec')==18
    
        protocolName{k,1} = 'CC_Capacitance';
        CCcount = CCcount + 1;                                              % count capacitance
        
    elseif sum(h.protocolName(end-12:end-4) == 'vclampRNA')==9
    
        protocolName{k,1} = 'VClampRNA';                                    % count vclamp RNA
        VCcount = VCcount + 1; 
        
    elseif sum(h.protocolName(end-12:end-4) == 'noise10pA')==9 ...
            || sum(h.protocolName(end-12:end-4) == 'noise30pA')==9 ...
            || sum(h.protocolName(end-12:end-4) == 'noise50pA')==9 ...
            || sum(h.protocolName(end-12:end-4) == 'noise70pA')==9 ...
            || sum(h.protocolName(end-12:end-4) == 'noise90pA')==9 ...
            || sum(h.protocolName(end-13:end-4) == 'noise110pA')==10 ...
            || sum(h.protocolName(end-13:end-4) == 'noise130pA')==10 ...
            || sum(h.protocolName(end-13:end-4) == 'noise150pA')==10 ...
            || sum(h.protocolName(end-13:end-4) == 'noise170pA')==10 ...
            || sum(h.protocolName(end-13:end-4) == 'noise190pA')==10 ...
            || sum(h.protocolName(end-13:end-4) == 'noise210pA')==10
        
        protocolName{k,1} = 'noise';
        Noisecount = Noisecount + 1;                                        % count noise
        
    else
        
        protocolName{k,1} = 'protocol name could not be identified';
        NotAvailcount = NotAvailcount + 1;                                  % count where N/A
        
    end
end

% count to see if all files were analyzed for protocol info
fullCount = (LPcount + SPcount + GFcount + SSPcount + RampUpcount + ...
    NAcount + RampNewcount + Noisecount + NotAvailcount + CCcount + ...
    VCcount - 11) / length(files);                                          % ==1 if all files