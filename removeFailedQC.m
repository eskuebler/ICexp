%{
removeFailedQC
%}

fieldnames_var = fieldnames(IC);                                            % Getting the variable names to overwrite in them
fieldnames_var([78,90]) = [];

for n = 1:length(cellList)  
    if isnan(IC.resistance_hd(n,1)) || isnan(IC.rheobaseLP(n,1))            % If crucial features cannot be determined, all parameters are set to NaN  
        for var = 11:length(fieldnames_var)-2                               % Leave the first 7 and last 2 untouched, since they are still usefull
            IC.(fieldnames_var{var})(n,1) = NaN;
        end
%         IC.peakTimeLP{n,1} = NaN;
    end
 end
