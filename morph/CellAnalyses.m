%
clear; close all; clc
start_trees
%free parameters
dd = 10; %Difference in diameter (um)
title = 'please find your data folder';
dataFolder = uigetdir('C:\',title);
a = dir(dataFolder);
a = a(2:end);
fbar = waitbar(0,'Analyzed cells','Name','Progress','CreateCancelBtn' ,'setappdata(gcbf,''canceling'',1)');
for n=1:length(a)
    if getappdata(fbar,'canceling')
        break
    end
    waitbar(n/length(a),fbar)
    dirStr = [a(n,1).folder,'/',a(n,1).name];
    b = dir(dirStr);
    b = b(3:end);
    if length(b) == 22
        %neurolucida_tree ('Test .ASC','-s')
        trees = load_tree ([b(2,1).folder,'/',b(2,1).name],'-s');
        
        Analysis
        %export_fig('C:/Users/joshuapoulin/Documents/MATLAB/Trees/Output', '-png', '-transparent')
     
        close
        clear trees intree
    end
end
delete(fbar)