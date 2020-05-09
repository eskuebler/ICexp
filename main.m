%{
intracellular experiment pipeline
%}

% Matlab startup commands
clear; close all; clc;

preSP = 20; postSP = 75;
preLP = 600; postLP = 2000;
plotSweeps = 1;
missingCountLP = 1; missingCountSP = 1;
 
% enter folder where data resides (each cell needs a folder)
dataFolder{1} = 'C:\Users\EKuebler\OneDrive - The University of Western Ontario\NHP cell type database (JMT)\01 main';
dataFolder{2} = 'D:\cell data';

start_trees; clear trees

% generate structures for each cell (i.e., ephys and morph)
for m = 1:length(dataFolder)                                                % folders denoting dataset
    cellList = dir(dataFolder{m}); cellList = cellList(3:end);
    Norig = length(cellList);                                               % number of cells in the dataset
    for n = 1:Norig
        disp(cellList(n).name)
        fileList = dir([dataFolder{m},'\',cellList(n,1).name,'/']);         % could contain various protocols or morphology
        fileList = fileList(3:end);
        LPcount = 1; SPcount = 1;
        for k = 1:length(fileList)
            exten = fileList(k,1).name(end-2:end);
            if sum(exten == 'abf')==3
                getSweepsABF
            elseif sum(exten == 'nwb')==3
                getSweepsNWB
            elseif sum(exten == 'swc')==3
                getMorphSWC
            end
            clear exten
        end
        
        if exist('LP') && exist('SP') && exist('trees')
            if plotSweeps == 1
                figure('Position',[50 50 1200 300]); set(gcf,'color','w');
                if LP.fullStruct == 1
                    for swp = 1:size(LP.V,2)
                        subplot(1,3,1); hold on; plot(LP.V{1,swp},'k','linewidth',0.25); axis tight; xlabel('time'); ylabel('V')
                    end
                else
                    disp(' stim info missing LP')
                    missingListLP{missingCountLP} = cellList(n).name;
                    missingCountLP = missingCountLP + 1;
                end
                if SP.fullStruct == 1
                    for swp = 1:size(SP.V,2)
                        subplot(1,3,2); hold on; plot(SP.V{1,swp},'k','linewidth',0.25); axis tight; xlabel('time'); ylabel('V')
                    end
                else
                    disp(' stim info missing SP')
                    missingListSP{missingCountSP} = cellList(n).name;
                    missingCountSP = missingCountSP + 1;
                end
                subplot(1,3,3); hold on; grid on;
                if length(trees)>1
                    for i = 1:length(trees)
                        ind = find(trees{1,i}.D==0);
                        if ~isempty(ind)
                            vec = 1:length(trees{1,i}.D);
                            vec(ind) = [];
                            trees{1,i}.D(ind,1) = min(trees{1,i}.D(vec,1));
                        end
                        scatter3(trees{1,i}.X,trees{1,i}.Y,trees{1,i}.Z,trees{1,i}.D,'k'); view(-36,59); axis tight;
                    end
                else
                    scatter3(trees.X,trees.Y,trees.Z,trees.D,'k'); view(-36,59); axis tight;
                end
                export_fig(['D:\genpath\',cellList(n).name,' raw data'],'-pdf','-r100');
                close
            end
            save(['D:\genpath\',cellList(n,1).name,'.mat'],'LP','SP','trees'); clear LP SP trees
        elseif ~exist('LP') && exist('SP') && exist('trees')
            if plotSweeps == 1
                figure('Position',[50 50 800 300]); set(gcf,'color','w');
                if SP.fullStruct == 1                    
                    for swp = 1:size(SP.V,2)
                        subplot(1,2,2); hold on; plot(SP.V{1,swp},'k','linewidth',0.25); axis tight; xlabel('time'); ylabel('V')
                    end
                else
                    disp(' stim info missing SP')
                    missingListSP{missingCountSP} = cellList(n).name;
                    missingCountSP = missingCountSP + 1;
                end
                subplot(1,2,2); hold on; grid on;
                if length(trees)>1
                    for i = 1:length(trees)
                        ind = find(trees{1,i}.D==0);
                        if ~isempty(ind)
                            vec = 1:length(trees{1,i}.D);
                            vec(ind) = [];
                            trees{1,i}.D(ind,1) = min(trees{1,i}.D(vec,1));
                        end
                        scatter3(trees{1,i}.X,trees{1,i}.Y,trees{1,i}.Z,trees{1,i}.D,'k'); view(-36,59); axis tight;
                    end
                else
                    scatter3(trees.X,trees.Y,trees.Z,trees.D,'k'); view(-36,59); axis tight;
                end
                export_fig(['D:\genpath\',cellList(n).name,' raw data'],'-pdf','-r100');
                close
            end
            save(['D:\genpath\',cellList(n,1).name,'.mat'],'SP','trees'); clear SP trees
        elseif exist('LP') && ~exist('SP') && exist('trees')
            if plotSweeps == 1
                figure('Position',[50 50 800 300]); set(gcf,'color','w');
                if LP.fullStruct  == 1
                    for swp = 1:size(LP.V,2)
                        subplot(1,2,1); hold on; plot(LP.V{1,swp},'k','linewidth',0.25); axis tight; xlabel('time'); ylabel('V')
                    end
                else
                    disp(' stim info missing LP')
                    missingListLP{missingCountLP} = cellList(n).name;
                    missingCountLP = missingCountLP + 1;
                end
                subplot(1,2,2); hold on; grid on;
                if length(trees)>1
                    for i = 1:length(trees)
                        ind = find(trees{1,i}.D==0);
                        if ~isempty(ind)
                            vec = 1:length(trees{1,i}.D);
                            vec(ind) = [];
                            trees{1,i}.D(ind,1) = min(trees{1,i}.D(vec,1));
                        end
                        scatter3(trees{1,i}.X,trees{1,i}.Y,trees{1,i}.Z,trees{1,i}.D,'k'); view(-36,59); axis tight;
                    end
                else
                    scatter3(trees.X,trees.Y,trees.Z,trees.D,'k'); view(-36,59); axis tight;
                end
                export_fig(['D:\genpath\',cellList(n).name,' raw data'],'-pdf','-r100');
                close
            end
            save(['D:\genpath\',cellList(n,1).name,'.mat'],'LP','trees'); clear LP trees
        elseif exist('LP') && exist('SP') && ~exist('trees')
            if plotSweeps == 1
                figure('Position',[50 50 800 300]); set(gcf,'color','w');
                if LP.fullStruct == 1
                    for swp = 1:size(LP.V,2)
                        subplot(1,2,1); hold on; plot(LP.V{1,swp},'k','linewidth',0.25); axis tight; xlabel('time'); ylabel('V')
                    end
                else
                    disp(' stim info missing LP')
                    missingListLP{missingCountLP} = cellList(n).name;
                    missingCountLP = missingCountLP + 1;
                end
                if SP.fullStruct == 1
                    for swp = 1:size(SP.V,2)
                        subplot(1,2,2); hold on; plot(SP.V{1,swp},'k','linewidth',0.25); axis tight; xlabel('time'); ylabel('V')
                    end
                else
                    disp(' stim info missing SP')
                    missingListSP{missingCountSP} = cellList(n).name;
                    missingCountSP = missingCountSP + 1;
                end
                export_fig(['D:\genpath\',cellList(n).name,' raw data'],'-pdf','-r100');
                close
            end
            save(['D:\genpath\',cellList(n,1).name,'.mat'],'LP','SP'); clear LP SP
        elseif exist('LP') && ~exist('SP') && ~exist('trees')
            if plotSweeps == 1
                figure('Position',[50 50 400 300]); set(gcf,'color','w');
                if LP.fullStruct == 1
                    for swp = 1:size(LP.V,2)
                        hold on; plot(LP.V{1,swp},'k','linewidth',0.25); axis tight; xlabel('time'); ylabel('V')
                    end
                else
                    disp(' stim info missing LP')
                    missingListLP{missingCountLP} = cellList(n).name;
                    missingCountLP = missingCountLP + 1;
                end
                export_fig(['D:\genpath\',cellList(n).name,' raw data'],'-pdf','-r100');
                close
            end
            save(['D:\genpath\',cellList(n,1).name,'.mat'],'LP'); clear LP
        elseif exist('SP') && ~exist('LP') && ~exist('trees')
            if plotSweeps == 1
                figure('Position',[50 50 400 300]); set(gcf,'color','w');
                if SP.fullStruct == 1
                    for swp = 1:size(SP.V,2)
                        hold on; plot(SP.V{1,swp},'k','linewidth',0.25); axis tight; xlabel('time'); ylabel('V')
                    end
                else
                    disp(' stim info missing SP')
                    missingListSP{missingCountSP} = cellList(n).name;
                    missingCountSP = missingCountSP + 1;
                end
                export_fig(['D:\genpath\',cellList(n).name,' raw data'],'-pdf','-r100');
                close
            end
            save(['D:\genpath\',cellList(n,1).name,'.mat'],'SP'); clear SP
        end
        clear k fileList
    end
    clear n cellList Norig
end
