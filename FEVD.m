function [VarDec] = FEVD(Theta, B, Lambda, spec, T, h, State)

[Theta_SR_N, ~] = CalculateImpulseResponses(spec, Theta, B,  T, h);

Vec = ones(T(1,2)^2, 1 );
for i = 2: T(1,2)^2
    Vec(i) = Vec(i-1) + 1;
end
Vec = reshape ( reshape(Vec, T(1,2), T(1,2))', T(1,2)^2, 1 )';          % positions of IRs, example {1 4 7 2 5 8 . 6 9 }

L = diag(get_sigma(Lambda, T, State) );         % diagonal of Lambda 
L_dup = kron(sqrt(L), ones(T(1,2), 1));                              % weight of the Lambda

% MSE computation
MSE=zeros(T(1,2) , h);

for j = 1:h
    for i = 1: T(1,2)   % var
       Mse = 0;
       for k = 1:T(1,2)
            Mse = Mse +  ( ( sqrt(L(k)) * Theta_SR_N( Vec( (i-1)*T(1,2) + k ), j) ) ^2 ) ;
       end
       if j > 1
            MSE(i, j) = MSE(i, j-1) +  Mse;
       else
           MSE(i, j) =  Mse;
       end     
    end

end

% computation of each component
Theta_Sq(:, 1) = ( (L_dup .* Theta_SR_N(:, 1) ) .^ 2 ) ;
for j = 2:h
    for i = 1: T(1,2)^2   % var
            Theta_Sq(i, j) =  Theta_Sq(i, j-1) + ( (L_dup(i) * Theta_SR_N(i, j) )  .^2)   ;
     end
end


% Cmputation of the contribution of each shock to total variance
for  j = 1:h
    cnt = 1;
    for i = 0 : T(1, 2) : ( T(1, 2)^2 - 1 );
        for k = 1:T(1,2)
            VarDec(i+k, j) = Theta_Sq( Vec(i+k), j) / MSE(cnt, j);
        end
        cnt = cnt + 1;
    end
end


