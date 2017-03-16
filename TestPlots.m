figure
plot(squeeze(yMAatPeriodTbyShocks(2,1, :)) )
hold on;
plot(squeeze(yMAatPeriodTbyShocks(2,2, :)) )

figure
plot(data(7: end,2) )
hold on;
plot( data(7: end,2) - squeeze(yMAatPeriodTbyShocks(2,2, :)) )

figure
plot(yAccumulated(7: end,2) )
hold on;
plot( yAccumulatedMAatPeriodTbyShocks )

figure
plot(squeeze(yMAatPeriodTbyShocks(2,1, :)))
hold on;
plot(squeeze(yMAatPeriodTbyShocks(2,2, :)))

figure
plot(yAccumulated(7: end,2)  - yAccumulatedMAatPeriodTbyShocks(:,2))
hold on;
plot(yAccumulated(7: end,2) )

figure
plot(yMA(:,1)  - squeeze(yMAatPeriodTbyShocks(1,2,:)), 'k')
hold on;
plot(yMA(:,1) )

figure1 = figure('Name','Counterfactual', 'NumberTitle','off', 'PaperType','<custom>','PaperSize', [11.81 4.5]);
positions = 1: size(yAccumulatedMAatPeriodTbyShocks, 1);
plot(yAccumulated(7: end,2), 'k', 'LineWidth', 1.5)
hold on;
maximum = 2;
fill([0 positions size(yAccumulatedMAatPeriodTbyShocks, 1)], [maximum yAccumulated(7: end,2)' maximum], [0.9 0.9 0.9], 'EdgeColor', 'k')
hold on;
plot( yAccumulatedMAatPeriodTbyShocks(:,2), 'k--', 'LineWidth', 1.5)
fill([0 positions size(yAccumulatedMAatPeriodTbyShocks, 1)], [maximum yAccumulatedMAatPeriodTbyShocks(:,2)' maximum], [0.8 0.8 0.8], 'EdgeColor', 'none')
hold off;
set(gca, 'XTickLabel',{'1995','2000','2005','2010','2015'}, 'XTick',sort( 290: -60 : 1) );
set(gca, 'XLim', [1 302]);
legend('\pi_t^{e, l}','Conbtribution of news shocks', 'Target shocks in the historical decomposition' ,'Contribution of target shocks');

figure
%plot(squeeze(yMAatPeriodTbyShocks(2,2, :)))
hold on;
plot(yAccumulatedMA(:,2), 'k')
plot(squeeze(yAccumulatedMAatPeriodTbyShocks(:, 2)), 'k--')
plot(squeeze(yAccumulatedMAatPeriodTbyShocks(:, 3)), 'b--')
plot(squeeze(yAccumulatedMAatPeriodTbyShocks(:, 2)) + squeeze(yAccumulatedMAatPeriodTbyShocks(:, 3)), 'b')
