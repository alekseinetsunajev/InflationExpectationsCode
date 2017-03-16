% This function estimates the Impulse Responses and computes various
% bootstrap Confidence Intervals. Procedure follows Lütkepohl (2005) D.3,
% p.709.
function [Theta_SR, Theta_LR, Bootstrap, Bootstrap_LR,  bootstrap_Y_MA] = Bootstrap(spec, Theta, B, B_stdErr, KsiT, Lambda, T, y, Z, h, brep)

% To compute IRs with actual data:
if B(1,1) < 0
    B(:,1) = -1 * B(:,1);
end
if B(2,2) < 0
    B(:,2) = -1 * B(:,2);
end

[Theta_SR, Theta_LR] = CalculateImpulseResponses(spec, Theta, B,  T, h);

U = GetResiduals(T, y, Z, spec.lags, Theta);
C = Theta(1:T(1,2))' ;

%% The bootstrap loop
Bootstrap = zeros(T(1,2) ^ 2,  (h+1), brep);
for n=1:brep
    fprintf(1,'Current bootstrap replication ');
    fprintf(1,'%d \n', n);
     % Set initial values:
    y_star = zeros(T(1,1), T(1,2) );
    y_star(1: spec.lags, :) = y(1:spec.lags, :); % Pre-sample values    
    
    %% Fixed design wild bootstrap series 
    
    % Rademacher distribution:
    rade = 2*round(rand(T(1,1)-spec.lags, 1))-1;
    
    %  use the Rademacher distribution,
    for t=1:(T(1,1)-spec.lags)
        Ustar(t, :) = (U(t, :)' * rade(t,:) )';
    end
          
    % Generate wild bootstraped series,
    for t=spec.lags+1 : T(1,1)
        if spec.trend == 1
            for i=T(1,2):T(1,2):T(1,2)*spec.lags
                y_star(t, :)  = y_star(t, :) + (get_coefficient(Theta, T, i/T(1,2), spec) * y(t-i/T(1,2), :)')' ;
            end
            y_star(t, :) = y_star(t, :) + C + t * Theta(T(1,2)+1:2*T(1,2))'  + Ustar(t-spec.lags, :);
            
        else
            for i=T(1,2):T(1,2):T(1,2)*spec.lags
                y_star(t, :)  = y_star(t, :) + (get_coefficient(Theta, T, i/T(1,2), spec) * y(t-i/T(1,2), :)')' ;
            end
            y_star(t, :) = y_star(t, :) + C + Ustar(t-spec.lags, :);
        end
    end

    % Obtain parameters based on bootstrapped series    
    Y = y_star(spec.lags+1:end, :);
    A1 = 0;
    A2 = 0;
    for i = 1:T(1,1)-spec.lags;
        Z_star(i,1) = 1;            % constant
        if spec.trend == 1      % with trend
            Z_star(i,2) = i;
            for j=1:spec.lags
                Z_star(i,(j-1)*T(1,2)+2+1 : j*T(1,2)+1+1 )= y(spec.lags + i - j , :);  % lagged observations 
            end
        else
             for j=1:spec.lags
                Z_star(i,(j-1)*T(1,2)+2 : j*T(1,2)+1 )= y(spec.lags + i - j , :);  % lagged observations 
            end           
        end
        
        A1 = A1 + kron( Z_star(i,:)' * Z_star(i,:), eye(T(1,2)));         
        A2 = A2 + kron(Z_star(i,:)', eye(T(1,2)))*y_star(i+spec.lags,:)'; 

    end
    
    RND=1;
    val(1,1)=1;
    val(1,2)=1;
    
    Tm = sum(KsiT ) - KsiT(1,:);
   Ahat1_vec  = Theta;
   U_opt =  GetResiduals(T, y_star, Z_star, spec.lags, Ahat1_vec);

    for position = 1:spec.s
        Sigma(1 + T(1,2)*(position-1) : T(1,2) * position, : )= B * ...
                            Lambda(1 + T(1,2)*(position-1) : T(1,2) * position, :) * B' ;
    end
    L_full(:, 1) = reshape(Lambda,  T(1,2)^2*spec.s, 1 );

    while val(RND,2) > 0.0000000001
      
       [BL, val(RND+1,1), xi] = MinimizeN(KsiT, T, Tm, Lambda, B, spec, U_opt, Ahat1_vec, 1);

        % Store B 
        B = BL(:,1: T(1,2));
        
        % store percentage change
        val(RND+1,2) = abs( ( val(RND+1,1) - val(RND,1) ) / val(RND,1) );

        % store lambda
         for position = 2 : spec.s
             Lambda(1 + T(1,2)*(position-1) : T(1,2) * position, :) = diag( BL(:, T(1,2)+position-1 ) );
         end

        if B(1,1) < 0
                B(:,1) = -1 * B(:,1);
        end
        if B(2,2) < 0
            B(:,2) = -1 * B(:,2);
        end
        % Write the covariance matrices from the decomposition of B and Lambda
       for position = 1:spec.s
          Sigma (1 + T(1,2)*(position-1) : T(1,2) * position, :) = B * ...
                Lambda(1 + T(1,2)*(position-1) : T(1,2) * position, :) * B';      
       end
        % Estimate parameter vector Theta (Ahat1_vec)
        T2=0;                                   %auxiliary counter matrices
        T3=0;                                   %auxiliary counter matrices
        for position = 1:spec.s
            T2 = T2 + kron( ( Z_star' * diag(KsiT(2:T(1,1)-spec.lags+1,position)) * Z_star ), get_sigma(Sigma, T, position)^-1) ;
            T3 = T3 + kron( Z_star' * diag(KsiT(2:T(1,1)-spec.lags+1,position)), get_sigma(Sigma, T, position)^-1);
        end
        Ahat1_vec = T2^-1 * T3 * reshape(Y', (T(1,1)-spec.lags) * T(1,2),1) ;
        U_opt = GetResiduals(T, y_star, Z_star, spec.lags, Ahat1_vec);        
        RND = RND + 1;
    end
        
    B_full(:, n) = reshape(B, T(1,2)^2, 1);
    L_full(:, n+1) = reshape(Lambda,  T(1,2)^2*spec.s, 1 );    
      
    [Theta1, Theta1_LR] = CalculateImpulseResponses(spec, Ahat1_vec, B,  T, h);
    % Store the bootstrapped IR values in matrix
    Bootstrap(:,:, n) =Theta1 ;
    Bootstrap_LR(:,:, n) =  Theta1_LR;
    
    % calculate FEVD(Theta, B, Lambda...)
    FEVD1(:,:, n) = FEVD(Ahat1_vec, B,  Lambda, spec, T, 100, 1);
    FEVD2(:,:, n) = FEVD(Ahat1_vec, B,  Lambda, spec, T, 100, 2);
    
    estimationBootstrap.Theta = Ahat1_vec;
    estimationBootstrap.B = B;
    bootstrap_Y_MA(:,:, n) = HistDecSimulation(estimationBootstrap, T, y, Z, spec);    
end

%% Compute Confidence Intervals
q = 0.975;
CIH2 = quantile(Bootstrap_LR,q,3);
CIL2 = quantile(Bootstrap_LR,1-q,3);

% rescale the impulse responses s.t. the shocks have a unit impact
Theta = Theta_LR;
Theta([1, 3],:) = Theta_SR([1, 3],:);
CIH2([1, 3], :) = quantile(Bootstrap([1, 3],:,:),q,3);
CIL2([1, 3], :) = quantile(Bootstrap([1, 3],:,:),1-q,3);

cnt = 0;
%pos = [2 1 4 3];        % position
pos = [1 2 3 4];        % position
for i = 1:T(1,2):T(1,2)^2
    CIL2(i: i + T(1,2) - 1, :) = CIL2(i : i + T(1,2) - 1, :) ./ Theta_SR(pos(i + cnt),1);
    CIH2(i : i + T(1,2) - 1, :) = CIH2(i : i + T(1,2) - 1, :) ./ Theta_SR(pos(i + cnt) ,1);
    Theta(i:i + T(1,2) - 1, :) = Theta(i:i + T(1,2) - 1, :) ./ Theta_SR(pos(i + cnt),1);
    cnt = cnt + 1;
end

%% Plot of Impulse Responses with Confidence Intervals
% PICTURE FOR INFLATION EXPECTATIONS
figure('Name','Impulse Responses','NumberTitle','off')
periods=0:h;
Vec = reshape(reshape((1:T(1,2)^2), T(1,2), T(1,2))', T(1,2)^2, 1)';
for p = 1 : T(1,2)^2
                subplot(T(1,2), T(1,2), p);
                plot(periods, Theta(Vec(p),:),'k',  'LineWidth', 1.5);
                hold on;
                plot(periods, CIH2(Vec(p),:),'k--',  'LineWidth', 1.5);
                plot(periods, CIL2(Vec(p),:),'k--' ,  'LineWidth', 1.5);
                plot(periods, zeros(1, h+1), ':' );
                hold off;
                xlim(gca, [0 h]);
end



