function [Z, A1, A2] = CalculateMatrixOfLaggedValues (y, lags, trend)
A1 = 0;
A2 = 0;
T = size(y);
Z = zeros(T(1,1)- lags, lags *T(1,2) + 1 );

for i = 1:T(1,1) - lags;
    Z(i,1) = 1;                                         % constant
    
    if trend == 1
        % with trend
        Z(i,2) = i;
        for j=1: lags
            Z(i,(j-1)*T(1,2)+2+1: j*T(1,2)+1+1 )= y(lags+i-j,:);  % lagged observations 
        end
    else
        % without trend
        for j=1: lags
            Z(i,(j-1)*T(1,2)+2: j*T(1,2)+1 )= y(lags + i-j,:);  % lagged observations 
        end            
    end
    A1 = A1 + kron( Z(i,:)' * Z(i,:), eye(T(1,2)));
    A2 = A2 + kron(Z(i,:)', eye(T(1,2)))*y(i + lags,:)'; 
end
