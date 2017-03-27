% Fumction that specifies and minimizes the lihelihood
% function where cov. matrices have structural decomposition on B and Lambda

function[BL, val] = MinimizeN(KsiT, T, Tm, Lambda, B, spec, u, Theta, ROUND)

options = optimset('Algorithm', 'interior-point', 'Display', 'off', 'MaxFunEvals', 20000, ...
    'MaxIter', 20000, 'Hessian','bfgs', 'DerivativeCheck','on','Diagnostics','off','GradObj','off','LargeScale','off', 'UseParallel', 'always');

% Form the matrix x of parameters to optimize wrt
parametersToOptimise = B;
for position = 1:spec.s-1
    parametersToOptimise(:,T(1,2)+position) = diag(Lambda(position*T(1,2)+1:T(1,2)*(position+1), :));
end

% lb is lower bound for constrained optimization, lbound returns the
%  bound for different number of states and variables
lb = lbound(T, spec);

% Restrictions  
restriction = [0 1 0 0];
% selection matrix for the sorting of Lambda elements
selectionMatrix = zeros(T(1,2)-1, T(1,2));
for k = 1 : (T(1,2)-1)
    selectionMatrix(k, k:k+1) = [1 -1];
end

if spec.s == 2
  
    W = eye(T(1,2),T(1,2));
    for w=1:spec.lags
        W = W - GetCoefficient(Theta(:,ROUND), T, w);
    end

    ineq = @(x)[ selectionMatrix * x(:, T(1,2)+1); ];
    %  Restrictions for long run 
    ceq = @(x) [
      restriction * reshape(W^-1*x(1:T(1,2), 1: T(1,2)), T(1,2)^2 , 1)
     ];
    
     if spec.BQrestrict == 1
        NonLinconstr = @(x)deal([] ,ceq(x));
     else
        NonLinconstr = @(x)deal( ineq(x), []);
     end
     
elseif spec.s == 3
          
    W = eye(T(1,2),T(1,2));
    for w=1:spec.lags
        W = W - GetCoefficient(Theta(:,ROUND), T, w);
    end

	ineq = @(x)[];
    % Restrictions for long run  
    ceq = @(x) [ 
    restriction * reshape(W^-1*x(1:T(1,2), 1: T(1,2)), T(1,2)^2 , 1)
    ];
    
    if spec.BQrestrict == 1
        NonLinconstr = @(x)deal( [ineq(x)],ceq(x));
    else
        NonLinconstr = @(x)deal( ineq(x), []);
    end
 
elseif spec.s == 4
    W = eye(T(1,2),T(1,2));
    for w=1:spec.lags
        W = W - GetCoefficient(Theta(:,ROUND), T, w);
    end
    
     ineq = @(x)[
       ];
    
    ceq = @(x) [ 
    restriction * reshape(W^-1*x(1:T(1,2), 1: T(1,2)), T(1,2)^2 , 1)
    ];
    
    if spec.BQrestrict == 1
        NonLinconstr = @(x)deal( [] ,ceq(x));
    else
        NonLinconstr =  @(x)deal( ineq(x), []);
    end
    
else
    disp('Number of states exceeds the one that is implemented (4) ');
    return
end
%% Optimization part
[BL, val] = fmincon( @(parametersToOptimise)LogLike(parametersToOptimise, KsiT, Tm, spec, u, T), parametersToOptimise , [], [] , [] , [] , lb , [] , NonLinconstr, options);

