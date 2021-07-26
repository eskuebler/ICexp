function [output,isis] = median_isi(LP)
isis = [];

for k = 1:length(LP.stats)
   if isfield(LP.stats{k, 1},'ISI')
       isis = [isis, LP.stats{k, 1}.ISI];
   end
end

output = 1000/(nanmedian(isis));
end