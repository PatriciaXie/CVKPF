classdef MOTScore
    %MOTScore 对跟踪进行打分
    %   此处显示详细说明
    
    properties
        rsTrackers
        gtTrackers
        MOTA
        MOTP
        Recall
        Precision
        IDSW
    end
    
    methods
        function ms = MOTScore(gtTrajectories, rsTrackers)
            %MOTScore 构造此类的实例
            % 1 真值
            nTrajectory = length(gtTrajectories);
            trackers = [];
            id = 0;
            for i = 1:1:nTrajectory
                pieces = getPieces(gtTrajectories{i}(:, 9));
                for k = 1:1:size(pieces, 1)
                    tracker.start = pieces(k, 1);
                    tracker.end = pieces(k, 2);
                    id = id + 1;
                    tracker.id = id;
                    tracker.states = [];
                    for j = tracker.start:1:tracker.end
                            tracker.states = [tracker.states, gtTrajectories{i}(j, :)'];
                    end
                    trackers = [trackers; tracker];
                    clear tracker
                end
            end
            ms.gtTrackers = trackers;
            clear trackers
            
            ms.rsTrackers = [];
            % 2 对跟踪进行编号
            for i = 1:1:length(rsTrackers)
                tracker.id = i;
                tracker.start = rsTrackers(i).start;
                tracker.end = rsTrackers(i).end;
                tracker.states = rsTrackers(i).states;
                ms.rsTrackers = [ms.rsTrackers; tracker];
            end
        end
        
        function outputFile(ms, gtFilename, rsFilename)
            %outputFile 将tracker结构的真值和跟踪输出为
            %   此处显示详细说明
            % gt
            fid=fopen(gtFilename,'w');
            for i = 1:1:length(ms.gtTrackers)
                tracker = ms.gtTrackers(i);
                for j = 1:1:(tracker.end - tracker.start + 1)
                    iframe = tracker.start + j - 1;
                    id = tracker.id;
                    x = tracker.states(1, j);
                    y = tracker.states(2, j);
                    z = tracker.states(3, j);
                    theta = tracker.states(4, j);
                    phi = tracker.states(5, j);
                    fprintf(fid,' %d, %d, %.5f, %.5f, %.5f, %.5f, %.5f\n', iframe, id, x, y, z, theta, phi);
                end
            end
            % rs
            fclose(fid);
            fid=fopen(rsFilename,'w');
            for i = 1:1:length(ms.rsTrackers)
                tracker = ms.rsTrackers(i);
                for j = 1:1:(tracker.end - tracker.start + 1)
                    iframe = tracker.start + j - 1;
                    id = tracker.id;
                    x = tracker.states(1, j);
                    y = tracker.states(2, j);
                    z = tracker.states(3, j);
                    theta = tracker.states(4, j);
                    phi = tracker.states(5, j);
                    fprintf(fid,' %d, %d, %.5f, %.5f, %.5f, %.5f, %.5f\n', iframe, id, x, y, z, theta, phi);
                end
            end
            fclose(fid);
        end
        
        function ms = getScore(ms, box)
            % gt
            gtMat = [];
            for i = 1:1:length(ms.gtTrackers)
                tracker = ms.gtTrackers(i);
                for j = 1:1:(tracker.end - tracker.start + 1)
                    iframe = tracker.start + j - 1;
                    id = tracker.id;
                    x = tracker.states(1, j);
                    y = tracker.states(2, j);
                    z = tracker.states(3, j);
                    theta = tracker.states(4, j);
                    phi = tracker.states(5, j);
                    gtMat = [gtMat; iframe, id, x, y, z, theta, phi];
                end
            end
            
            % rs
            rsMat = [];
            for i = 1:1:length(ms.rsTrackers)
                tracker = ms.rsTrackers(i);
                for j = 1:1:(tracker.end - tracker.start + 1)
                    iframe = tracker.start + j - 1;
                    id = tracker.id;
                    x = tracker.states(1, j);
                    y = tracker.states(2, j);
                    z = tracker.states(3, j);
                    theta = tracker.states(4, j);
                    phi = tracker.states(5, j);
                    rsMat = [rsMat; iframe, id, x, y, z, theta, phi];
                end
            end
            
            threshold = 0.003;
            [ms.Recall, ms.Precision, ms.IDSW, ms.MOTA, ms.MOTP] = CLEAR_MOT(gtMat, rsMat, threshold);
            
            if ~isempty(box)
                figure(1)
                hold on
                xlabel('X');
                ylabel('Z');
                zlabel('Y');
                axis equal
                grid on
                view(-5, -10);
                plot3(1000 * box.Xb, 1000 * box.Yb, 1000 * box.Zb, 'k');
%                 for iFly = 1:1:length(ms.rsTrackers)
%                     X =  1000 * ms.rsTrackers(iFly).states(1, :);
%                     Y =  1000 * ms.rsTrackers(iFly).states(2, :);
%                     Z =  1000 * ms.rsTrackers(iFly).states(3, :);
%                     plot3(X, Z, Y, '-', 'LineWidth', 2);
%                     %     plot3(X(1), Y(1), Z(1), 'o');
%                 end
                %                 set(gca,'ZDir','reverse');
                %                 title('Track Result');
                %                 hold off
                
                %                 figure(2), clf
                %                 hold on
                %                 xlabel('X');
                %                 ylabel('Z');
                %                 zlabel('Y');
                %                 axis equal
                %                 grid on
                %                 view(-5, -10);
                %                 plot3(box.Xb, box.Yb, box.Zb, 'k');
                for iFly = 1:1:length(ms.gtTrackers)
                    X =  1000 * ms.gtTrackers(iFly).states(1, :);
                    Y =  1000 * ms.gtTrackers(iFly).states(2, :);
                    Z =  1000 * ms.gtTrackers(iFly).states(3, :);
                    plot3(X, Z, Y); %, 'k.');
                    %     plot3(X(1), Y(1), Z(1), 'o');
                end
                set(gca,'ZDir','reverse');
                info = sprintf('MOTA: %.2f %%  MOTP: %.2f mm\n', ms.MOTA, ms.MOTP);
                title(info);
                hold off
%                 set(gcf,'color','none');
%                 set(gca,'color','none'); 
            end
            
        end
    end
end

