function[Y_bar_ma] = HistDecSimulation(Estimation, T, y, Z, spec)
% obtain historical decomposition

% calculate struct. shocks
U = GetResiduals(T, y, Z, spec.lags, Estimation.Theta);

for ii = 1:T(1, 1) - spec.lags
    Epsilon(ii, :) = (Estimation.B^-1 * U(ii, :)');
end

% convert to VAR(1)
A = zeros(T(1, 2)* (spec.lags) ); % Define the A matrix
for i = 1:spec.lags
      A(1:T(1,2), T(1,2)*(i-1)+1 : T(1,2)*i ) = get_coefficient(Estimation.Theta, T, i, spec);
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
     y_bar = y(spec.lags + 1, :)';

    for jj = 1 : ii        
        y_bar = y_bar + J*A^(jj-1)*J' * Estimation.B * counterfactualEpsilon(ii - jj +1, :)';
    end
    y_bar_ma(ii, :) =  y_bar' ;
   y_bar_ma(ii, 2) = roundn( y_bar_ma(ii, 2), -1);    
end

% accumulate data
y_A(1, :) = [3.6 4];
y_A = [y_A; y];
for  ii = 2:T(1, 1) +1
    y_A(ii, :) = y_A( ii-1, :) + y_A( ii, :); 
end
% accumulate artificial data out of shocks. NB! rounding is very important!
Y_bar_ma = [y_A(spec.lags+1, :) ; y_bar_ma];
for  ii = 2:T(1, 1) +1- spec.lags
    Y_bar_ma(ii, :)  = Y_bar_ma(ii-1, :) + Y_bar_ma(ii, :);
    %Y_bar_ma(ii, :)  = roundn( Y_bar_ma(ii, :), -1);
end
y_A = y_A(2:end, :);
Y_bar_ma = Y_bar_ma(2:end, :);

