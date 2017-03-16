%%
%Function produces standard errors, returns - StErrors is the array of
%standard errors that have been estimated by method specified below in the variable Method,
%StErr_all - array of standard errors where first column presents the
%finite difference standard errors, second column - the outer product
%matrix standard errors and third column - White standard errors (see below)

function[StErrors, StErr_all, xi] = Param_Stderr(T, Estimation, spec, y, Z)
%% Choose the likelihood function
Tm = sum(Estimation.KsiT ) - Estimation.KsiT(1,:);
j=1;

if spec.s == 2
    param = [reshape(Estimation.B, T(1,2)^2,1)', diag(Estimation.Lambda(T(1,2)+1 : T(1,2)*2,:))', Estimation.Theta', reshape(Estimation.P(1:spec.s-1, :), numel(Estimation.P(1:spec.s-1, :)), 1 )'];
%    param = [reshape(Estimation.B, T(1,2)^2,1)', diag(Estimation.Lambda(T(1,2)+1 : T(1,2)*2,:))', Estimation.Theta'];

    Elements = numel(param);
    Elem_Theta = numel(Estimation.Theta);
    Elem_P = numel(Estimation.P(1:spec.s-1, :));
    LogLVec = @(x)(T(1,1)-spec.lags)*log( abs( det( reshape(x(1:T(1,2)^2, 1), T(1,2), T(1,2)) ) ) ) + 0.5*trace( reshape(x(1:T(1,2)^2, 1), T(1, 2), T(1, 2))'^-1 * reshape(x(1:T(1,2)^2, 1), T(1, 2),T(1, 2))^-1 ...
    * (repmat(Estimation.KsiT(2:T(1,1) - spec.lags + 1, 1),1, T(1,2) ) .* GetResiduals(T, y, Z, spec.lags, x(Elements - Elem_Theta - Elem_P + 1 : Elements - Elem_P) ) )' * GetResiduals(T, y, Z, spec.lags, x(Elements - Elem_Theta - Elem_P + 1 : Elements - Elem_P)) ) + ...
            Tm(1,2)/2 * log( det( diag( x(T(1,2)^2 +1 : T(1,2)^2 + T(1,2), 1) ) ) ) + 0.5*trace( reshape(x(1:T(1,2)^2, 1), T(1,2), T(1,2))'^-1 * diag( x(T(1,2)^2+1 : T(1,2)^2 + T(1,2), 1) )^-1 * reshape(x(1:T(1,2)^2, 1), T(1, 2), T(1, 2))^-1 ...
        * (repmat(Estimation.KsiT(2:T(1,1) - spec.lags + 1, 2),1, T(1,2) ) .* GetResiduals(T, y, Z, spec.lags, x(Elements - Elem_Theta - Elem_P + 1 : Elements - Elem_P)) )' * GetResiduals(T, y, Z, spec.lags, x(Elements - Elem_Theta - Elem_P + 1 : Elements - Elem_P)) ) ;

elseif spec.s == 3
    param = [reshape(Estimation.B, T(1,2)^2,1)', diag(Estimation.Lambda(T(1,2)+1 : T(1,2)*2,:))', diag(Estimation.Lambda(T(1,2)*2+1:T(1,2)*3,:))', Estimation.Theta', reshape(Estimation.P(1:spec.s-1, :), numel(Estimation.P(1:spec.s-1, :)), 1 )'];
    Elements = numel(param);
    Elem_Theta = numel(Estimation.Theta);
    Elem_P = numel(Estimation.P(1:spec.s-1, :));
    LogLVec = @(x)(T(1,1)-spec.lags)*log( abs( det( reshape(x(1:T(1,2)^2, 1), T(1,2), T(1,2)) ) ) ) + 0.5*trace( reshape(x(1:T(1,2)^2, 1), T(1, 2), T(1, 2))'^-1 * reshape(x(1:T(1,2)^2, 1), T(1, 2),T(1, 2))^-1 ...
        * (repmat(Estimation.KsiT(2:T(1,1) - spec.lags + 1, 1),1, T(1,2) ) .* GetResiduals(T, y, Z, spec.lags, x(Elements - Elem_Theta - Elem_P + 1 : Elements - Elem_P) ) )' * GetResiduals(T, y, Z, spec.lags, x(Elements - Elem_Theta - Elem_P + 1 : Elements - Elem_P) ) )+ ...
        Tm(1,2)/2 * log( det( diag( x(T(1,2)^2 +1 : T(1,2)^2 + T(1,2), 1) ) ) ) + 0.5*trace( reshape(x(1:T(1,2)^2, 1), T(1,2), T(1,2))'^-1 * diag( x(T(1,2)^2 +1 : T(1,2)^2 + T(1,2), 1) )^-1 * reshape(x(1:T(1,2)^2, 1), T(1, 2), T(1, 2))^-1 ...
        * (repmat(Estimation.KsiT(2:T(1,1) - spec.lags + 1, 2),1, T(1,2) ) .* GetResiduals(T, y, Z, spec.lags, x(Elements - Elem_Theta - Elem_P + 1 : Elements - Elem_P)) )' * GetResiduals(T, y, Z, spec.lags, x(Elements - Elem_Theta - Elem_P + 1 : Elements - Elem_P) ) ) + ...
        Tm(1,3)/2 * log( det( diag( x( T(1,2)^2 + T(1,2) + 1:T(1,2)^2 + 2*T(1,2), 1) ) ) ) + 0.5*trace( reshape(x(1:T(1,2)^2, 1), T(1,2), T(1,2))'^-1 * diag( x( T(1,2)^2 + T(1,2) + 1:T(1,2)^2 + 2*T(1,2), 1) )^-1 * reshape(x(1:T(1,2)^2, 1), T(1, 2),T(1, 2))^-1 ... 
        * (repmat(Estimation.KsiT(2:T(1,1) - spec.lags + 1, 3),1, T(1,2) ) .* GetResiduals(T, y, Z, spec.lags, x(Elements - Elem_Theta - Elem_P + 1 : Elements - Elem_P)) )' * GetResiduals(T, y, Z, spec.lags, x(Elements - Elem_Theta - Elem_P + 1 : Elements - Elem_P) ) ) ;
    
elseif spec.s == 4
    param = [reshape(Estimation.B, T(1,2)^2,1)', diag(Estimation.Lambda(T(1,2)+1 : T(1,2)*2,:))', diag(Estimation.Lambda(T(1,2)*2+1:T(1,2)*3,:))', diag(Estimation.Lambda(T(1,2)*3+1:T(1,2)*4,:))', Estimation.Theta', reshape(Estimation.P(1:spec.s-1, :), numel(Estimation.P(1:spec.s-1, :)), 1 )' ];
    Elements = numel(param);
    Elem_Theta = numel(Estimation.Theta);
    Elem_P = numel(Estimation.P(1:spec.s-1, :));
    LogLVec = @(x)(T(1,1)-spec.lags)*log( abs( det( reshape(x(1:T(1,2)^2, 1), T(1,2), T(1,2)) ) ) ) + 0.5*trace( reshape(x(1:T(1,2)^2, 1), T(1, 2), T(1, 2))'^-1 * reshape(x(1:T(1,2)^2, 1), T(1, 2),T(1, 2))^-1 ...
        * (repmat(Estimation.KsiT(2:T(1,1) - spec.lags + 1, 1),1, T(1,2) ) .* GetResiduals(T, y, Z, spec.lags, x(Elements - Elem_Theta - Elem_P + 1 : Elements - Elem_P) ) )' * GetResiduals(T, y, Z, spec.lags, x(Elements - Elem_Theta - Elem_P + 1 : Elements - Elem_P)) ) + ...
        Tm(1,2)/2 * log( det( diag( x(T(1,2)^2 +1 : T(1,2)^2 + T(1,2), 1) ) ) ) + 0.5*trace( reshape(x(1:T(1,2)^2, 1), T(1,2), T(1,2))'^-1 * diag( x(T(1,2)^2 +1 : T(1,2)^2 + T(1,2), 1) )^-1 * reshape(x(1:T(1,2)^2, 1), T(1, 2), T(1, 2))^-1 ...
        * (repmat(Estimation.KsiT(2:T(1,1) - spec.lags + 1, 2),1, T(1,2) ) .* GetResiduals(T, y, Z, spec.lags, x(Elements - Elem_Theta - Elem_P + 1 : Elements - Elem_P)) )' * GetResiduals(T, y, Z, spec.lags, x(Elements - Elem_Theta - Elem_P + 1 : Elements - Elem_P)) ) + ...
        Tm(1,3)/2 * log( det( diag( x( T(1,2)^2 + T(1,2) + 1:T(1,2)^2 + 2*T(1,2), 1) ) ) ) + 0.5*trace( reshape(x(1:T(1,2)^2, 1), T(1,2), T(1,2))'^-1 * diag( x( T(1,2)^2 + T(1,2) + 1:T(1,2)^2 + 2*T(1,2), 1) )^-1 * reshape(x(1:T(1,2)^2, 1), T(1, 2),T(1, 2))^-1 ... 
        * (repmat(Estimation.KsiT(2:T(1,1) - spec.lags + 1, 3),1, T(1,2) ) .* GetResiduals(T, y, Z, spec.lags, x(Elements - Elem_Theta - Elem_P + 1 : Elements - Elem_P)) )' * GetResiduals(T, y, Z, spec.lags, x(Elements - Elem_Theta - Elem_P + 1 : Elements - Elem_P)) ) + ... 
        Tm(1,4)/2 * log( det( diag( x( T(1,2)^2 + 2*T(1,2) + 1:T(1,2)^2 + 3*T(1,2), 1) ) ) ) + 0.5*trace( reshape(x(1:T(1,2)^2, 1), T(1,2), T(1,2))'^-1 * diag( x( T(1,2)^2 + 2*T(1,2) + 1:T(1,2)^2 + 3*T(1,2), 1) )^-1 * reshape(x(1:T(1,2)^2, 1), T(1, 2),T(1, 2))^-1 ... 
        * (repmat(Estimation.KsiT(2:T(1,1) - spec.lags + 1, 4),1, T(1,2) ) .* GetResiduals(T, y, Z, spec.lags, x(Elements - Elem_Theta - Elem_P + 1 : Elements - Elem_P)) )' * GetResiduals(T, y, Z, spec.lags, x(Elements - Elem_Theta - Elem_P + 1 : Elements - Elem_P)) ) ;


end
param = param';
myDelta=1e-4*abs(param)+0.00000001;

%%
% First derivative calculation as in Hamilton page 143 and Krolzig 120-121.
xi=1;
s=zeros(T(1,1)-spec.lags,numel(param));

% calculate the initial likelihood function for each t.
for t=1:T(1,1)-spec.lags
    logLikVec1(t,1) = log( ( Estimation.P * Estimation.Ksi(t,:)' )' * Estimation.Eta(t,:)');
end

for i=1:numel(param)
    m=param;
    m(i)=param(i)+myDelta(i);    
    
    B = reshape(m(1:T(1,2)^2, 1), T(1,2), T(1,2)) ;
    % Form matrix Sigma of state-dependent covariance matrices
    if spec.s == 2
        Sigma(1:T(1,2) , :) = B*B' ;
        Sigma(T(1,2)+1 : T(1,2)*2, :) = B * diag( m(T(1,2)^2 +1 : T(1,2)^2 + T(1,2), 1) ) * B' ;
    elseif spec.s == 3
        Sigma(1:T(1,2) , :) = B*B' ;
        Sigma(T(1,2)+1 : T(1,2)*2, :) = B * diag( m(T(1,2)^2 +1 : T(1,2)^2 + T(1,2), 1) ) * B' ;
        Sigma(2*T(1,2)+1 : T(1,2)*3, :) = B * diag( m(T(1,2)^2 + T(1,2) + 1 : T(1,2)^2 + 2*T(1,2), 1) ) * B' ;
    elseif spec.s == 4
        Sigma(1:T(1,2) , :) = B*B' ;
        Sigma(T(1,2)+1 : T(1,2)*2, :) = B * diag( m(T(1,2)^2 +1 : T(1,2)^2 + T(1,2), 1) ) * B' ;
        Sigma(2*T(1,2)+1 : T(1,2)*3, :) = B * diag( m(T(1,2)^2 + T(1,2) + 1 : T(1,2)^2 + 2*T(1,2), 1) ) * B' ;
        Sigma(3*T(1,2)+1 : T(1,2)*4, :) = B * diag( m(T(1,2)^2 + 2*T(1,2) + 1 : T(1,2)^2 + 3*T(1,2), 1) ) * B' ;
    end
    
	P = reshape(m( numel(m) - numel(Estimation.P(1:spec.s-1, :)) + 1 : numel(m)), spec.s-1 , spec.s);
	
    for k=1:spec.s-1
        for q=1:spec.s
            if P(k,q)>1
                P(k,q) = 0.999999;
            elseif P(k,q)<=0
                P(k,q) = 0.0000001;
            end
        end    
     end
    
    P_last = ones(1, spec.s);
	for k=1:spec.s-1
        P_last = P_last - P(k,:);
    end
	P(spec.s, :) = P_last;

    u = GetResiduals(T, y, Z, spec.lags, m(Elements - Elem_Theta - numel(Estimation.P(1:spec.s-1, :)) + 1 : Elements - numel(Estimation.P(1:spec.s-1, :))) );
    Ksi2 = Estimation.Ksi;

    % update conditional density function Eta for the changes in parameter
    % values
    for j = 1: T(1,1)-spec.lags
        for position=1:spec.s
            Eta2(j,position) = (2*pi)^(-T(1,2)/2)*det(get_sigma(Sigma, T, position) )^-0.5 * exp(-0.5*u(j,:)...
                *get_sigma(Sigma, T, position)^-1*u(j,:)' );
        end
    end    

    % update state probabilities for the changes in parameter values
    for q=2:T(1,1)-spec.lags + 1
        Ksi2(q,1:spec.s) = ( Eta2(q-1,:)' .* (P*Ksi2(q-1,1:spec.s)') /...
            (spec.const'*(Eta2(q-1,:)'.*(P*Ksi2(q-1,1:spec.s)'))) )'; 
    end
    
    % calculate the approximated derivative for each t and each parameter
    for t=1:T(1,1)-spec.lags
        logLikVec2(t,i) = log( (P * Ksi2(t,:)' )' *  Eta2(t,:)');
        s(t,i)=(logLikVec2(t,i)-logLikVec1(t,1))/(myDelta(i));
    end
   
end

% form the outer product matrix as in Hamilton  eq 5.8.4
sum_s_Matrix=zeros(numel(param));
for idx=1:T(1,1)-spec.lags
    s_Matrix=s(idx,:)'*s(idx,:);
    sum_s_Matrix=sum_s_Matrix+s_Matrix;
end
OP_Matrix=sum_s_Matrix/T(1,1);
StErr_all(:, 2) = sqrt(diag(inv(OP_Matrix*T(1,1))));


%%
%Produce results
% Using second partial derivatives, as eq 5.8.3 in Hamilton p 143 
%StErr_all(:, 1) = sqrt(diag(inv(H*T(1,1))));
% Using outer product Matrix, as eq 5.8.4 and below in Hamilton p 143
StErr_all(:, 2) = sqrt(diag(inv(OP_Matrix*T(1,1))));
%Using White's Covariance Matrix, as eq 5.8.7 in Hamilton p 145
StErr_all(:, 3) = ones(numel(param),1); %sqrt(diag(1/T(1,1)*((H*(OP_Matrix)^-1*H)^-1)));

method = 2;         % outer product matrix is the default method
switch method
    case 1  
        StErrors = StErr_all(:, 1);
    case 2  
        StErrors = StErr_all(:, 2);
    case 3  
        StErrors = StErr_all(:, 3);
end
