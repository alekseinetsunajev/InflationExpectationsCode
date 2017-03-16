function[] = BarGraph(dataToPlot)
figure
axes;
axis([1 size(dataToPlot, 1)+1 -0.5 4.65]);

color1 = [0.5 0.5 0.5];
color2 = [0.75 0.75 0.75];
% color1 = [255/255 255/255 0/255];
% color2 = [135/255 208/255 235/255];
hold on;

for ii = 1:size(dataToPlot, 1)
    if (dataToPlot(ii,1) > 0 && dataToPlot(ii,2) < 0 )
       h1 = rectangle('Position', [ii 0 1 dataToPlot(ii,1) ], 'FaceColor', color1);
       h2 = rectangle('Position', [ii dataToPlot(ii,2) 1 abs(dataToPlot(ii,2)) ], 'FaceColor', color2);
    end
    if (dataToPlot(ii,1) < 0 && dataToPlot(ii,2) > 0 )
        rectangle('Position', [ii dataToPlot(ii,1) 1 abs(dataToPlot(ii,1)) ], 'FaceColor', color1);
        rectangle('Position', [ii 0 1 dataToPlot(ii,2) ], 'FaceColor', color2);
    end
    if (dataToPlot(ii,1) < 0 && dataToPlot(ii,2) < 0 )
        rectangle('Position', [ii dataToPlot(ii,1) 1 abs(dataToPlot(ii,1)) ], 'FaceColor', color1);
        rectangle('Position', [ii dataToPlot(ii,2) + dataToPlot(ii,1)  1 abs(dataToPlot(ii,2)) ], 'FaceColor', color2);
    end
    if (dataToPlot(ii,1) > 0 && dataToPlot(ii,2) > 0 )
        h1 = rectangle('Position', [ii 0 1 abs(dataToPlot(ii,1)) ], 'FaceColor', color1);
        h2 = rectangle('Position', [ii dataToPlot(ii,1)  1 abs(dataToPlot(ii,2)) ], 'FaceColor', color2);
        p1 = plot(nan, nan, 's', 'markeredgecolor',get(h1,'edgecolor'), 'markerfacecolor',get(h1,'facecolor'));
        p2 = plot(nan, nan, 's', 'markeredgecolor',get(h2,'edgecolor'), 'markerfacecolor',get(h2,'facecolor'));
    end
    if (dataToPlot(ii,1) == 0 && dataToPlot(ii,2) < 0 )
       h2 = rectangle('Position', [ii dataToPlot(ii,2) 1 abs(dataToPlot(ii,2)) ], 'FaceColor', color2);
    end    
    if (dataToPlot(ii,1) == 0 && dataToPlot(ii,2) > 0 )
        rectangle('Position', [ii 0 1 dataToPlot(ii,2) ], 'FaceColor', color2);
    end    
end
%pd1 = plot(dataToPlot(:, 3), 'b',  'LineWidth',1.5);
pd2 = plot(dataToPlot(:, 4), 'k',  'LineWidth',1.5);
legend([p1 p2 pd2], 'News shocks', 'Target shocks', '\pi_t^{e,l}');





