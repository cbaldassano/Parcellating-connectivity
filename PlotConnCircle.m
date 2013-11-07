function PlotConnCircle(subjConn, subjLabels)


    s=1;
    segmentEndPts = zeros(12,1);
    gap = 15;
    smallgap=0;
    for k=1:12
        if (k<=2)
            nGaps = 0;
        elseif (k<=4)
            nGaps = 1;
        elseif (k<=6)
            nGaps = 2;
        elseif (k<=8)
            nGaps = 2.5;
        elseif (k<=10)
            nGaps = 3.5;
        elseif (k<=12);
            nGaps = 4;
        end
        segmentEndPts(k,1) = (find(subjLabels{s}==k,1,'first')-1+nGaps*gap+(k-1)*smallgap)/(length(subjLabels{s}==k)+5*gap+11*smallgap)*2*pi;
        segmentEndPts(k,2) = (find(subjLabels{s}==k,1,'last')-1+nGaps*gap+(k-1)*smallgap)/(length(subjLabels{s}==k)+5*gap+11*smallgap)*2*pi;
    end
    

%     breaks = linspace(0,2*pi,13);
%     segmentEndPts = [breaks(1:end-1)' breaks(2:end)'];

    radius = 5;
    centerpull = 0.05;
    numPts = 40;
    pullStrength = -4*centerpull/(numPts-1)*(((0:(numPts-1)).^2)/(numPts-1) - (0:(numPts-1)));
    pullStrength = (repmat(pullStrength',1,2)).^(0.5);
    
    figure('Color',[1 1 1],'Position',[ 2583         735         757         698]); hold on;

    colors = [0 186 99; 0 101 255]./255;
    for i=1:2:size(segmentEndPts,1)
        segmentAngles = linspace(segmentEndPts(i,1),segmentEndPts(i+1,2),100);
        h=plot(radius*cos(segmentAngles),radius*sin(segmentAngles),'-'); hold on;
        if (i<=4)
            set(h,'Color',colors(ceil(i/2),:));
        else
            set(h,'Color','k');
        end
        set(h,'LineWidth',3);
    end
    
    
    for s=1%:20
        for k = 4:-1:1%1:4
            segIndStart = find(subjLabels{s}==k,1,'first');
            segIndEnd = find(subjLabels{s}==k,1,'last');
            for i=segIndStart:segIndEnd
                for k2 = 5:12
                    segIndStart2 = find(subjLabels{s}==k2,1,'first');
                    segIndEnd2 = find(subjLabels{s}==k2,1,'last');
                    for j=segIndStart2:segIndEnd2
                        if (subjConn{s}(i,j)>0)
                            angleOut = segmentEndPts(k,1) + (segmentEndPts(k,2)-segmentEndPts(k,1))*(i-segIndStart)/(segIndEnd-segIndStart);
                            angleIn = segmentEndPts(k2,1) + (segmentEndPts(k2,2)-segmentEndPts(k2,1))*(j-segIndStart2)/(segIndEnd2-segIndStart2);
                            curve = [linspace(radius*cos(angleOut),radius*cos(angleIn),numPts)' linspace(radius*sin(angleOut),radius*sin(angleIn),numPts)'];
                            curve = (1-pullStrength).*curve;

                            h= plot(curve(:,1),curve(:,2));
                            set(h,'Color',colors(ceil(k/2),:));
                        end
                    end
                end
            end
        end
    end
    hold off;
    axis square; axis off;
            
        
end