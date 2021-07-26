function [sp] = estimatePeak(LP,int4Peak2,k)

%{
there is an issue here whereby several time points are the same voltage as
the detected 'peak' - to deal with this there is an if statement using
equalV to find last time point where the voltage was equal
%}

% estimate peak of AP
% if length(LP.putSpTimes2)==1
% 	for i = 1:length(int4Peak2)
% 		[peak(i), peakTime(i)]= max(LP.V{1,k}(int4Peak2{i}(1)-1:int4Peak2{i}(end)-1));
%         
% 	end
% elseif length(LP.putSpTimes2)>1
    for i = 1:length(int4Peak2)
        [peak(i), peakTime(i)]= max(LP.V{1,k}(int4Peak2{i}(1)-1:int4Peak2{i}(end)));
        equalV = LP.V{1,k}(int4Peak2{i}(1)-1:int4Peak2{i}(end))==peak(i);
        if sum(equalV)>1
            peakTime(i) = find(equalV,1,'last');
        end
    end
% end
peakTime = LP.putSpTimes2 + peakTime - 2;

sp.peak = peak;
sp.peakTime = peakTime;