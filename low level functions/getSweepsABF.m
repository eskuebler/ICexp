%{
getSweepsABF
%}

listToCorrect = {'M05_SM_A1_C02','M05_SM_A1_C03','M05_SM_A1_C06',...
    'M05_SM_A1_C07','M05_SM_A1_C09','M05_SM_A1_C13','M05_SM_A1_C15',...
    'M06_SM_A1_C01','M06_SM_A1_C06','M06_SM_A1_C07','M06_SM_A1_C08',...
    'M06_SM_A1_C09','M06_SM_A1_C12','M06_SM_A1_C14'};                   % stim times off

[d,~,h] = abfload([fileList(k).folder,'/',fileList(k).name], ...
    'sweeps','a','channels','a');
d = squeeze(d);
LP.fullStruct = 1; SP.fullStruct = 1;
if sum(h.protocolName(end-22:end-4) == ...
        'Monkey_1000 ms step')==19 ...
    || length(h.protocolName)>27 && sum(h.protocolName(end-26:end-4) == ...
    'Monkey_1000 ms steplong')==23

%% meagan diff protocol names for long pulse
    parametersABF_LP
elseif sum(h.protocolName(end-19:end-4) == 'Monkey_3 ms step')==16 ...
    || sum(h.protocolName(end-20:end-4) == 'Monkey_3 ms step2')==17 ...
    || sum(h.protocolName(end-28:end-4) == 'Monkey_3 ms step_addsweep')==25
    
    parametersABF_SP
end
if LPcount == 1
    LP.fullStruct = 0;
end
if SPcount == 1
    SP.fullStruct = 0;
end