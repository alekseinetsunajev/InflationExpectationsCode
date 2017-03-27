function [varianceDecomposition] = FEVD(Theta, B, Lambda, spec, T, h, regime, levelsOfSecondVariable)

[shortRunResponses, accumulatedResponses] = CalculateImpulseResponses(spec, Theta, B,  T, h);

if levelsOfSecondVariable == 1
    shortRunResponses = [shortRunResponses(1, :) ; accumulatedResponses(2, :); shortRunResponses(3, :);  accumulatedResponses(4, :)] ;
end

positonOfResponses = ones(T(1,2)^2, 1 );
for i = 2: T(1,2)^2
    positonOfResponses(i) = positonOfResponses(i-1) + 1;
end
positonOfResponses = reshape ( reshape(positonOfResponses, T(1,2), T(1,2))', T(1,2)^2, 1 )';          % positions of IRs, example {1 4 7 2 5 8 . 6 9 }

shockVariance = diag(GetSigmaForRegime(Lambda, T, regime) );         % diagonal of Lambda 
weightOfShockVariance = kron(sqrt(shockVariance), ones(T(1,2), 1));                              % weight of the Lambda

% MSE computation
MSE=zeros(T(1,2) , h);

for j = 1:h
    for i = 1: T(1,2)   % var
       Mse = 0;
       for k = 1:T(1,2)
            Mse = Mse +  ( ( sqrt(shockVariance(k)) * shortRunResponses( positonOfResponses( (i-1)*T(1,2) + k ), j) ) ^2 ) ;
       end
       if j > 1
            MSE(i, j) = MSE(i, j-1) +  Mse;
       else
           MSE(i, j) =  Mse;
       end     
    end

end

% computation of each component
contributionOfEachCompionent(:, 1) = ( (weightOfShockVariance .* shortRunResponses(:, 1) ) .^ 2 ) ;
for j = 2:h
    for i = 1: T(1,2)^2   % var
            contributionOfEachCompionent(i, j) =  contributionOfEachCompionent(i, j-1) + ( (weightOfShockVariance(i) * shortRunResponses(i, j) )  .^2)   ;
     end
end


% Cmputation of the contribution of each shock to total variance
for  j = 1:h
    cnt = 1;
    for i = 0 : T(1, 2) : ( T(1, 2)^2 - 1 );
        for k = 1:T(1,2)
            varianceDecomposition(i+k, j) = contributionOfEachCompionent( positonOfResponses(i+k), j) / MSE(cnt, j);
        end
        cnt = cnt + 1;
    end
end


