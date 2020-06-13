function output = median_isi(LP)
isis = [];
k = [];

for k = 1:length(LP.stats)
   if  isfield(LP.stats{k, 1},'spTimes') 
    if LP.stats{k, 1}.spTimes > 0
       isis = [isis, LP.stats{k, 1}.ISI];
    end
   end
end

output = 1000/(median(isis));

close