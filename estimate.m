%% Function that estimates the parameters of the model using the EM
% algorithm
function [ThetaLast, Sigma, KsiT, P, IterationLik, B, Lambda, Eta, Ksi] = estimate(A1, A2, T, y, Z, spec, Tolvalue, RandNumbers, LM)
%%  Set initial parameters
	RND_ALL = 1;
    
    P = 1/spec.s * spec.const * spec.const';
    
    % Get initial estimates of parameter vector and calculate initial residual matrix u
    Theta =  (A1)^-1 * A2 ;
    u = residuals2(T, y, Z, spec, Theta);           %function retrieves residuals

    %SET THE INITIAL PARAMETERS
    B =  (1/T(1,1)*(u'*u))^0.5 + RandNumbers;
    Ksi = (spec.const/spec.s)';

    % Form Lambda matrices for decomposition, the first Lambda_1 is identity
    % matrix
    for position = 1:spec.s
        Lambda(1 + T(1,2)*(position-1) : T(1,2) * position, :) = (LM^(position-1)) *eye(T(1,2));
    end

    % Form initial covariance matrices
    for position = 1:spec.s
        Sigma(1 + T(1,2)*(position-1) : T(1,2) * position, : )= B * ...
                            Lambda(1 + T(1,2)*(position-1) : T(1,2) * position, :) * B' ;
    end

    % Form initial likelihood values
    LogValue(1,1) = 0;
    LogValue(1,2) = 1;
    LogValue(1,3) = 1;

    % END OF INITIALIZATION STAGE


    % EM ALGORITHM RUNS UNTIL CONVERGENCE (while loop)
    % ROUND is a counter to accumulate number of iterations
    ROUND = 1;
    while abs(LogValue(ROUND, 3)) > Tolvalue        
%%      %EXPECTATION STEP

        %Calculate the conditional distribution function Eta   
        for k = 1:T(1,1)-spec.lags
            % function get_sigma retrieves covariance matrix from
            % Estimation.Sigma according to the position (state of the system)
            for position=1:spec.s
                Eta(k,position) = (2*pi)^(-T(1,2)/2)*det(get_sigma(Sigma, T, position) )^-0.5 * exp(-0.5*u(k,:)...
                            *get_sigma(Sigma, T, position)^-1*u(k,:)' );
            end 
        end

        % Filtering probabilities    
        for m=2:T(1,1)-spec.lags+1
            Ksi(m,1:spec.s) = ( Eta(m-1,:)' .* (P*Ksi(m-1,1:spec.s)') /...
                (spec.const'*(Eta(m-1,:)'.*(P*Ksi(m-1,1:spec.s)'))) )'; 
        end
        
        % avoid label switching
        Ksi(T(1,1)-spec.lags,:) = sort( Ksi(T(1,1)-spec.lags,:) , 2);
        
        % Smoothing probabilities    
        KsiTRevTime = Ksi(T(1,1)-spec.lags+1,:) ;                   %Inverse time!

        counter = T(1,1)-spec.lags+1;
        for n=1:T(1,1)-spec.lags
            counter = counter - 1;
            KsiTRevTime(n+1,1:spec.s) = ( (P' * (KsiTRevTime(n,1:spec.s)' ./ ...
            (P*Ksi(counter,1:spec.s)'))) .* Ksi(counter,1:spec.s)' )';
        end 

        KsiT = fliplr(rot90(rot90(KsiTRevTime)));       % Flip the matrix to obtain timing 0,1,2..T-1
        Tm = sum(KsiT ) - KsiT(1,:);                            % Calculate sum of probabilities

        for o=1:T(1,1)-spec.lags
            KsiT2(o,1:spec.s^2) = ... 
           (reshape(P', spec.s^2,1) .* kron ((KsiT(o+1,1:spec.s)'  ./ (P*Ksi(o,1:spec.s)')) , ...
                                                Ksi(o,1:spec.s)')  )' ;
        end

 %%     % Compute the value of likelihood function (computation is done for the previous
        % iteration) and the change in the likelihood 
        LogL = 0;
        for t=1:T(1,1)-spec.lags
            LogL = LogL + log( ( P * Ksi(t,:)' )' * Eta(t,:)') ;
        end
        LogValue(ROUND+1,1) = ROUND;
        LogValue(ROUND+1,2) = LogL;
        LogValue(ROUND+1,3) = ( LogValue(ROUND+1,2) - LogValue(ROUND,2) ) / LogValue(ROUND,2) ;

%%  MAXIMIZATION STEP

        %   Estimate the Markov matrix, the first inpup in the reshape() function is the
        %   vectorized Markov matrix
        P = reshape( sum(KsiT2)' ./ kron ( spec.const , kron(spec.const', eye(spec.s)) * sum(KsiT2)' ),...
             spec.s, spec.s)';

        % Minimize the likelihood function and retrieve the matricex B and
        % Lambdas, BL stores  B first, then it has diagonal entries of Lambdas using
        % function Minimize
        val = 0;
        val(1,1)=1;
        val(1,2)=1;
        RND = 1;
        while val(RND,2) > 0.000000001
            u = residuals2(T, y, Z, spec, Theta(:, RND_ALL));
            [BL, val(RND+1,1), ~] = MinimizeN(KsiT, T, Tm, Lambda, B, spec, u, Theta, RND_ALL);
            
            val(RND+1,2) = abs( ( val(RND+1,1) - val(RND,1) ) / val(RND,1) );

            % Store B
            B = BL(:,1: T(1,2));
            % Store Lambdas 
            for position = 2:spec.s
                Lambda(1 + T(1,2)*(position-1) : T(1,2) * position, :) = diag( BL(:, T(1,2)+position-1 ) );
            end

            % Write the covariance matrices from the decomposition of B and Lambda
            for position = 1:spec.s
                Sigma (1 + T(1,2)*(position-1) : T(1,2) * position, :) = B * ...
                    Lambda(1 + T(1,2)*(position-1) : T(1,2) * position, :) * B';      
            end

            % Estimate parameter vector Theta
            T2=0;                                   %auxiliary counter matrices
            T3=0;                                   %auxicounter matrices
            for position = 1:spec.s
                T2 = T2 + kron( ( Z' * diag(KsiT(2:T(1,1)-spec.lags+1,position)) * Z ), get_sigma(Sigma, T, position)^-1) ;
                T3 = T3 + kron( Z' * diag(KsiT(2:T(1,1)-spec.lags+1,position)), get_sigma(Sigma, T, position)^-1);
            end
            Theta(:, RND_ALL+1) = T2^-1 * T3 * reshape(y(spec.lags+1:T(1,1),:)', (T(1,1)-spec.lags) * T(1,2),1) ;
            RND = RND + 1;
            RND_ALL =RND_ALL + 1;
        end
        % update initial p
        Ksi(1,:) = KsiT(1,:);   

        fprintf(1,'Current round of iteration of EM is ');
        fprintf(1,'%d \n', ROUND);

        ROUND = ROUND + 1;
        
    end
       
ThetaLast = Theta(:, RND_ALL); 
IterationLik = [LogValue(ROUND,2), val(RND,1)];