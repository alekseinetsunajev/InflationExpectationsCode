function[] = Hist_dec_rnd(EstimationUR, T, y, Z, spec, CounterfactualDataQuantile)
% obtain historical decomposition

% calculate struct. shocks
U = residuals2(T, y, Z, spec, EstimationUR.Theta);

for ii = 1:T(1, 1) - spec.lags
    Epsilon(ii, :) = (EstimationUR.B^-1 * U(ii, :)');
end


% convert to VAR(1)
A = zeros(T(1, 2)* (spec.lags) ); % Define the A matrix
for i = 1:spec.lags
      A(1:T(1,2), T(1,2)*(i-1)+1 : T(1,2)*i ) = get_coefficient(EstimationUR.Theta, T, i, spec);
end
for m=1:spec.lags-1
    A(m*T(1,2)+1 : m*T(1,2)+T(1,2), T(1,2)*m-T(1,2)+1 : T(1,2) * m) = eye(T(1,2));
end
J = [eye(T(1,2)) zeros(T(1,2), T(1,2)*spec.lags - T(1,2) )  ];

% produce counterfactual series of errors
counterfactualEpsilon = zeros(T(1,1) - spec.lags, T(1,2));
counterfactualEpsilon(:, 1) = Epsilon(:, 1);

counterfactualEpsilon(:, 2) = Epsilon(:, 2);
counterfactualEpsilon(215:end, 1) = zeros(88, 1);

% MA representation   
for ii = 1:T(1,1) - spec.lags
    %y_bar = zeros(T(1,2),1);
     y_bar = y(spec.lags+1, :)';

    for jj = 1 : ii        
        y_bar = y_bar + J*A^(jj-1)*J' * EstimationUR.B * counterfactualEpsilon(ii - jj +1, :)';
    end

    y_bar_ma(ii, :) =  y_bar' ;
end

% accumulate data, obtain true series for the second variable
y_A(1, :) = [3.6 4];
y_A = [y_A; y];
for  ii = 2:T(1, 1) +1
    y_A(ii, :) = y_A( ii-1, :) + y_A( ii, :); 
end
% accumulate artificial data out of shocks. NB! rounding is very important!
Y_bar_ma = [y_A(spec.lags+1, :) ; y_bar_ma];
for  ii = 2:T(1, 1) +1- spec.lags
    Y_bar_ma(ii, :)  = Y_bar_ma(ii-1, :) + Y_bar_ma(ii, :);
    Y_bar_ma(ii, :)  = roundn( Y_bar_ma(ii, :), -1);
end
y_A = y_A(2:end, :);
Y_bar_ma = Y_bar_ma(2:end, :);

% plot the results
p = 2;

S(:,1) = y_A(spec.lags + 1 : T(1,1), p);
S(:,2) = 3*ones(T(1,1) - spec.lags , 1);
S(:,3) = Y_bar_ma(:, p);
S(:,4) = y_A(spec.lags + 1 : T(1,1), p)- Y_bar_ma(:, p);
S(:,5) = zeros(T(1,1) - spec.lags , 1);
S(:, 6) = CounterfactualDataQuantile.high(:, p);
S(:, 7) = CounterfactualDataQuantile.low(:, p);

figure1 = figure('Name','Counterfactual', 'NumberTitle','off', 'PaperType','<custom>','PaperSize',[11.81 4.5]);
[AX,H1,H2] = plotyy(1:1:88, S(215:end,[1 2 3 6 7]), 1:1:88, S(215:end,[4 5 ]) , 'plot');

set(get(AX(1),'Ylabel'),'String','Scale for \pi_t^{e, l} ', 'Interpreter', 'tex', 'color', [0 0 0], 'fontsize', 12) 
set(get(AX(2),'Ylabel'),'String','Scale for \{\pi_t^{e, l} - \pi_t^{e, l}counterfactual\}', 'Interpreter', 'tex', 'color', [0 0 0], 'fontsize', 12)
% change the X axes
axis(AX(1),[1 88 2 5]);               % axes for first set
axis(AX(2),[1 88 -1 0.3]);         % axes for second set

% assign the lables to X axes
set(AX(1), 'XTickLabel',{'2009','2010','2011','2012','2013','2014','2015'},...
    'XTick',[7 19 31 43 55 67 79 ], 'fontsize', 12, 'ycolor', [0 0 0], ...
    'YTickLabel',{'2.4', '2.6', '2.8', '3', '3.2', '3.4', '3.6'}, 'YTick', [2.4 2.6 2.8 3 3.2 3.4 3.6]);
set(AX(2), 'YTickLabel',{'-0.4', '-0.3', '-0.2', '-0.1', '0', '0.1'},'YTick',[-0.4 -0.3 -0.2 -0.1 0 0.1 ], 'XTickLabel',{},'XTick',[], 'fontsize', 12, 'ycolor', [0 0 0]);

% assign line styles and line width
set(H1(1), 'DisplayName','\pi_t^{e, l}', 'LineWidth',1.5,'Color',[0 0 0])
set(H1(3), 'DisplayName','\pi_t^{e, l} counterfactual', 'LineWidth',1.5,'Color',[0 0 0], 'LineStyle','--')
set(H1(2), 'DisplayName','', 'LineWidth',0.5,'Color',[0 0 0], 'LineStyle',':')
set(H1(4), 'DisplayName','95% bands for counterfactual', 'LineWidth',0.5,'Color',[0 0 0], 'LineStyle','-.')
set(H1(5), 'DisplayName','', 'LineWidth',0.5,'Color',[0 0 0], 'LineStyle','-.')

set(H2(1), 'DisplayName','\pi_t^{e, l} - \pi_t^{e, l}counterfactual', 'LineWidth',1.5,'Color',[0 0 0], 'LineStyle',':')
set(H2(2), 'DisplayName','', 'LineWidth',0.5,'Color',[0 0 0], 'LineStyle',':')

% % make the legend
 legend([H1(1); H1(3); H1(4); H2(1)]);
