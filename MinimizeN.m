% Fumction that specifies and minimizes the lihelihood
% function where cov. matrices have structural decomposition on B and Lambda

function[BL, val, xi] = MinimizeN(KsiT, T, Tm, Lambda, B, spec, u, Theta, ROUND)

options = optimset('Algorithm', 'interior-point', 'Display', 'off', 'MaxFunEvals', 20000, ...
    'MaxIter', 20000, 'Hessian','bfgs', ...
'DerivativeCheck','on','Diagnostics','off','GradObj','off','LargeScale','off', 'UseParallel', 'always');

% Form the matrix x of parameters to optimize wrt
x = B;
for position = 1:spec.s-1
    x(:,T(1,2)+position) = diag(Lambda(position*T(1,2)+1:T(1,2)*(position+1), :));
end

xV =  reshape(x, T(1,2)^2 + (spec.s-1)*T(1,2), 1);

% lb is lower bound for constrained optimization, lbound returns the
%  bound for different number of states and variables
lb = lbound(T, spec);

% Restrictions  
restriction = [0 1 0 0];




% selection matrix for the sorting of Lambda elements
select_matr = zeros(T(1,2)-1, T(1,2));
for k = 1 : (T(1,2)-1)
    select_matr(k, k:k+1) = [1 -1];
end

if spec.s == 2
  
    W = eye(T(1,2),T(1,2));
    for w=1:spec.lags
        W = W - get_coefficient(Theta(:,ROUND), T, w, spec);
    end

    ineq = @(x)[ select_matr * x(:, T(1,2)+1); ];

   ceq = @(x) [
   %  Restrictions for long run 
%     restriction * reshape( x(1:T(1,2), 1:T(1,2)), T(1,2)^2 , 1)   
      restriction * reshape(W^-1*x(1:T(1,2), 1: T(1,2)), T(1,2)^2 , 1)
     ];
    
    if spec.BQrestrict == 1
        NonLinconstr = @(x)deal([] ,ceq(x));
    else
        NonLinconstr = @(x)deal( ineq(x), []);
    end
    xi = 1;
    
elseif spec.s == 3
          
    W = eye(T(1,2),T(1,2));
    for w=1:spec.lags
        W = W - get_coefficient(Theta(:,ROUND), T, w, spec);
    end

	ineq = @(x)[];
    
    ceq = @(x) [ 
% Restrictions for long run
%     restriction * reshape( x(1:T(1,2), 1:T(1,2)), T(1,2)^2 , 1)
    restriction * reshape(W^-1*x(1:T(1,2), 1: T(1,2)), T(1,2)^2 , 1)
    ];
    
    if spec.BQrestrict == 1
        NonLinconstr = @(x)deal( [ineq(x)],ceq(x));
    else
        NonLinconstr = @(x)deal( ineq(x), []);
    end
 
xi = 1;
elseif spec.s == 4
    xi = 1;
    W = eye(T(1,2),T(1,2));
    for w=1:spec.lags
        W = W - get_coefficient(Theta(:,ROUND), T, w, spec);
    end
    
     ineq = @(x)[
       ];
    
    ceq = @(x) [
     ];
    
    if spec.BQrestrict == 1
        NonLinconstr = @(x)deal( [] ,ceq(x));
    else
        NonLinconstr =  @(x)deal( ineq(x), []);
    end
    
else
    'Number of states exceeds the one that is implemented (4) '
    return
end
%% Optimization part
[BL, val] = fmincon( @(x)LogLike(x, KsiT, Tm, spec, u, T), x , [], [] , [] , [] , lb , [] , NonLinconstr, options);
