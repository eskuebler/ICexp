    if length(class(trees))==6 & class(trees) == 'struct'
            intree = trees(1,1);
            options = '-s';
            sholl = sholl_tree   (intree, dd, options);
            s.Sholl{n,1} = sholl;
            [stats] = stats_tree   (intree, [], [], options);
            s.Stats{n,1} = stats;
            dstat = stats.dstats.BO{1, 1};
            s.Dstat{n,1} = dstat;
            cvol = cvol_tree    (intree, '-s');
            s.Cvol{n,1} = cvol;
            direct = direction_tree    (intree, '-s');
            s.Direct{n,1} = direct;
            eucl = eucl_tree    (intree, [], '-s',);
            s.Eucl{n,1} = eucl;
            len = len_tree (intree, '-s');
            Pvec = Pvec_tree (intree, len, options);
            Tort = Pvec./eucl;
            s.medianTort(n,1) = nanmedian(Tort);
          
          
            % SURF(X,Y,Z,C,...)
            surf = surf_tree    (intree, '-s');
            s.Surf{n,1} = surf;
            T = T_tree (intree);
            asym = asym_tree (intree, T, '');
            s.meanasym(n,1) = nanmean (asym);
            close all
            subplot(4,4,1:8)
            scatter3(intree.X,intree.Y, intree.Z, 4)
            subplot(4,4,9)
            plot(sholl)
            xlabel('distance from soma')
            ylabel('number of intersections')
            subplot(4,4,12)
            plot(direct)
            xlabel('direction')
            ylabel('amplitude')
            subplot(4,4,13)
            histogram(eucl)
            xlabel('euclidean distance')
            ylabel('count')
            subplot(4,4,14)
            histogram(Tort)
            xlabel('path length/euclidean distance')
            ylabel('count')
            subplot(4,4,15)
            histogram(cvol)
            xlabel('continuous volumes')
            ylabel('count')
            subplot(4,4,16)
            plot(surf)
            xlabel('surface values')
            ylabel('diameter (um^2)')
            clear sholl stats dstat angleB cvol len direct eucl surf
    end