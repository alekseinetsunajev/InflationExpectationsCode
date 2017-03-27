function[] = HistoricalDecomposition(Estimation, T, data, Z, spec)
% obtain historical decomposition

% calculate struct. shocks
residuals = GetResiduals(T, data, Z, spec.lags, Estimation.Theta);

for ii = 1:T(1, 1) - spec.lags
    structuralShocks(ii, :) = (Estimation.B^-1 * residuals(ii, :)');
end

% convert to VAR(1)
A = zeros(T(1, 2)* (spec.lags) ); % Define the A matrix
for i = 1:spec.lags
      A(1:T(1,2), T(1,2)*(i-1)+1 : T(1,2)*i ) = GetCoefficient(Estimation.Theta, T, i);
end
for m=1:spec.lags-1
    A(m*T(1,2)+1 : m*T(1,2)+T(1,2), T(1,2)*m-T(1,2)+1 : T(1,2) * m) = eye(T(1,2));
end
J = [eye(T(1,2)) zeros(T(1,2), T(1,2)*spec.lags - T(1,2) )  ];

% produce counterfactual series of errors
counterfactualStructuralShocks = zeros(T(1,1) - spec.lags, T(1,2));
counterfactualStructuralShocks(:, 1) = structuralShocks(:, 1);
counterfactualStructuralShocks(:, 2) = structuralShocks(:, 2);
%counterfactualStructuralShocks(215:end, 1) = zeros(88, 1);

% MA representation
initialY =  data(spec.lags, :)';
for ii = 2:spec.lags
	initialY((ii-1)*T(1,2) +1 :  (ii-1)*T(1,2) + T(1,2) , :) = data(spec.lags - ii + 1, :)';
end

mean = (eye(T(1,2) * spec.lags) - A)^-1 * [Estimation.Theta(1:T(1,2), 1 ); zeros(T(1,2) * (spec.lags -1 ), 1) ];
mean = mean(1:T(1,2));
yMAatPeriodTbyShocks  = zeros(T(1,2), T(1,2), T(1,1) - spec.lags);

for t = 1:T(1,1) - spec.lags
    initialA = A^t;
    yMAatPeriodT = initialA(1:T(1,2), :) * initialY + mean;
    constantTerm = yMAatPeriodT;
    yMAatPeriodTbyShocks(1,1, t)  = constantTerm(1,1);
    yMAatPeriodTbyShocks(2,2, t)  = constantTerm(2,1);
    for jj = 1 : t        
        yMAatPeriodT = yMAatPeriodT + J*A^(jj-1)*J' * Estimation.B * counterfactualStructuralShocks(t - jj +1, :)';
        
        matrixMA = J*A^(jj-1)*J' * Estimation.B;
        for variable = 1:T(1,2)
            for shock = 1:T(1,2)
                yMAatPeriodTbyShocks(variable,shock, t)  = yMAatPeriodTbyShocks(variable, shock, t) + matrixMA(variable, shock) * counterfactualStructuralShocks(t - jj +1, shock)';
            end
        end
    end
    yMA(t, :) = roundn(yMAatPeriodT', -1);
    yMATest(t, :) = sum(yMAatPeriodTbyShocks(:,:, t), 2 )';
end

firsDiffData = [squeeze(yMAatPeriodTbyShocks(2,1, :)) squeeze(yMAatPeriodTbyShocks(2,2, :))];
compareData = [data(spec.lags + 1:end, :) yMA yMATest];

% accumulate data, obtain true series for the second variable
yAccumulated(1, :) = [3.6 4];
yAccumulated = [yAccumulated; data];
for  ii = 2:T(1, 1) +1
    yAccumulated(ii, :) = yAccumulated( ii-1, :) + yAccumulated( ii, :); 
end

yAccumulated = yAccumulated(2:end, :);
startAccumulationAt = 5;
yAccumulatedMA = [yAccumulated(spec.lags + 1 : startAccumulationAt + spec.lags  - 1, :) ; yMA(startAccumulationAt : end,:)];
yAccumulatedMAatPeriodTbyShocks = yAccumulated(spec.lags + 1 : spec.lags  + startAccumulationAt - 1, :);
yAccumulatedMAatPeriodTbyShocks(startAccumulationAt-1, 3)  = yMAatPeriodTbyShocks(2,1, startAccumulationAt-1);

for  ii = startAccumulationAt : T(1, 1) - spec.lags
    yAccumulatedMA(ii, :)  = yAccumulatedMA(ii-1, :) + yAccumulatedMA(ii, :);
    yAccumulatedMAatPeriodTbyShocks(ii, 2)  = yAccumulatedMAatPeriodTbyShocks(ii-1, 2) + yMAatPeriodTbyShocks(2,2, ii) ;
    yAccumulatedMAatPeriodTbyShocks(ii, 2) = roundn(yAccumulatedMAatPeriodTbyShocks(ii, 2), -1);
    yAccumulatedMAatPeriodTbyShocks(ii, 3)  = yAccumulatedMAatPeriodTbyShocks(ii-1, 3) + yMAatPeriodTbyShocks(2,1, ii) ;   
    yAccumulatedMAatPeriodTbyShocks(ii, 3)  = roundn(yAccumulatedMAatPeriodTbyShocks(ii, 3), -1);
end

% yAccumulatedMAatPeriodTbyShocks(1, 3)  = 0.0;
% for  ii = 2 : T(1, 1) - spec.lags
%     yAccumulatedMAatPeriodTbyShocks(ii, 2)  = yAccumulatedMAatPeriodTbyShocks(ii-1, 2) + yMAatPeriodTbyShocks(2,2, ii) ;
%     yAccumulatedMAatPeriodTbyShocks(ii, 3)  = yAccumulatedMAatPeriodTbyShocks(ii-1, 3) + yMAatPeriodTbyShocks(2,1, ii) ;    
% end

compareAccumulatedData = [yAccumulated(spec.lags + 1:end, 2) yAccumulatedMA(:, 2) yAccumulated(spec.lags + 1:end, 2) - yAccumulatedMA(:, 2) ];
compareAccumulatedDataByShocks = [yAccumulatedMA(:, 2) yAccumulatedMAatPeriodTbyShocks(:, 2) yAccumulatedMAatPeriodTbyShocks(:, 3) yAccumulatedMAatPeriodTbyShocks(:, 2) + yAccumulatedMAatPeriodTbyShocks(:, 3)  ];
compareAccumulatedDataByShocks(:,5) = compareAccumulatedDataByShocks(:,1) - compareAccumulatedDataByShocks(:,4);

% plot the results
BarGraph(compareAccumulatedDataByShocks(:, [3 2 4 1]))
p = 2;
start = 0;
figure1 = figure('Name','Counterfactual', 'NumberTitle','off', 'PaperType','<custom>','PaperSize',[11.81 4.5]);
subplot(1,2,1, 'fontsize', 12)
plot( data(7+start:end, 1), 'k', 'LineWidth',1);
hold on;
plot(yMA(start + 1:end, 1), 'k:', 'LineWidth',1.5);

subplot(1,2,2, 'fontsize', 12)
plot( yAccumulated(7+start:end, p), 'k', 'LineWidth',1);
hold on;
plot(yAccumulatedMA(start + 1:end, p), 'k:', 'LineWidth',1.5);

% figure1 = figure('Name','Historical Decomposition', 'NumberTitle','off', 'PaperType','<custom>','PaperSize', [11.81 4.5]);
% positions = 1: size(yAccumulatedMAatPeriodTbyShocks, 1);
% plot(yAccumulated(7: end,2), 'k', 'LineWidth', 1.5)
% hold on;
% maximum = 2;
% fill([0 positions size(yAccumulatedMAatPeriodTbyShocks, 1)], [maximum yAccumulated(7: end,2)' maximum], [0.9 0.9 0.9], 'EdgeColor', 'k')
% hold on;
% plot( yAccumulatedMAatPeriodTbyShocks(:,2), 'k--', 'LineWidth', 1.5)
% fill([0 positions size(yAccumulatedMAatPeriodTbyShocks, 1)], [maximum yAccumulatedMAatPeriodTbyShocks(:,2)' maximum], [0.8 0.8 0.8], 'EdgeColor', 'none')
% hold off;
% set(gca, 'XTickLabel',{'1995','2000','2005','2010','2015'}, 'XTick',sort( 290: -60 : 1) );
% set(gca, 'XLim', [1 302]);
% legend('\pi_t^{e, l}','Conbtribution of news shocks', 'Target shocks in the historical decomposition' ,'Contribution of target shocks');

