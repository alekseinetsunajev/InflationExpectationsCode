%% Function that estimates the parameters of the model using the EM
% algorithm
function [thetaFinal, Sigma, shoothedRegimes, P, likelihoodFinal, B, Lambda, likelihoodForRegimes, filteredProbabilities] = estimate(theta, T, y, Z, spec, Tolvalue, RandNumbers, LM)
%%  Set initial parameters
	totalRounds = 1;
    
    % Get initial estimates of parameter vector and calculate initial residual matrix u
    P = 1/spec.s * spec.const * spec.const';
    u = GetResiduals(T, y, Z, spec.lags, theta);           %function retrieves residuals

    %SET THE INITIAL PARAMETERS
    B =  (1/T(1,1)*(u'*u))^0.5 + RandNumbers;
    filteredProbabilities = (spec.const/spec.s)';

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
    logValue(1,1) = 0;
    logValue(1,2) = 1;
    logValue(1,3) = 1;

    % END OF INITIALIZATION STAGE
    % EM ALGORITHM RUNS UNTIL CONVERGENCE (while loop)
    % ROUND is a counter to accumulate number of iterations
    outerRound = 1;
    while abs(logValue(outerRound, 3)) > Tolvalue        
%%      %EXPECTATION STEP

        %Calculate the conditional distribution function Eta   
        for k = 1:T(1,1)-spec.lags
            % function GetSigmaForRegime retrieves covariance matrix from
            % Estimation.Sigma according to the position (state of the system)
            for position=1:spec.s
                likelihoodForRegimes(k,position) = (2*pi)^(-T(1,2)/2)*det(GetSigmaForRegime(Sigma, T, position) )^-0.5 * exp(-0.5*u(k,:)...
                            * GetSigmaForRegime(Sigma, T, position)^-1*u(k,:)' );
            end 
        end

        % Filtering probabilities    
        for m=2:T(1,1)-spec.lags+1
            filteredProbabilities(m,1:spec.s) = ( likelihoodForRegimes(m-1,:)' .* (P*filteredProbabilities(m-1,1:spec.s)') /...
                (spec.const'*(likelihoodForRegimes(m-1,:)'.*(P*filteredProbabilities(m-1,1:spec.s)'))) )'; 
        end
        
        % avoid label switching
        filteredProbabilities(T(1,1)-spec.lags,:) = sort( filteredProbabilities(T(1,1)-spec.lags,:) , 2);
        
        % Smoothing probabilities    
        reversedTimeFilteredProbs = filteredProbabilities(T(1,1)-spec.lags+1,:) ;                   %Inverse time!

        counter = T(1,1)-spec.lags+1;
        for n=1:T(1,1)-spec.lags
            counter = counter - 1;
            reversedTimeFilteredProbs(n+1,1:spec.s) = ( (P' * (reversedTimeFilteredProbs(n,1:spec.s)' ./ ...
            (P*filteredProbabilities(counter,1:spec.s)'))) .* filteredProbabilities(counter,1:spec.s)' )';
        end 

        shoothedRegimes = fliplr(rot90(rot90(reversedTimeFilteredProbs)));       % Flip the matrix to obtain timing 0,1,2..T-1
        Tm = sum(shoothedRegimes ) - shoothedRegimes(1,:);                            % Calculate sum of probabilities

        for o=1:T(1,1)-spec.lags
            KsiT2(o,1:spec.s^2) = ... 
           (reshape(P', spec.s^2,1) .* kron ((shoothedRegimes(o+1,1:spec.s)'  ./ (P*filteredProbabilities(o,1:spec.s)')) , ...
                                                filteredProbabilities(o,1:spec.s)')  )' ;
        end

 %%     % Compute the value of likelihood function (computation is done for the previous
        % iteration) and the change in the likelihood 
        logL = 0;
        for t=1:T(1,1)-spec.lags
            logL = logL + log( ( P * filteredProbabilities(t,:)' )' * likelihoodForRegimes(t,:)') ;
        end
        logValue(outerRound+1,1) = outerRound;
        logValue(outerRound+1,2) = logL;
        logValue(outerRound+1,3) = ( logValue(outerRound+1,2) - logValue(outerRound,2) ) / logValue(outerRound,2) ;

%%  MAXIMIZATION STEP

        %   Estimate the Markov matrix, the first inpup in the reshape() function is the
        %   vectorized Markov matrix
        P = reshape( sum(KsiT2)' ./ kron ( spec.const , kron(spec.const', eye(spec.s)) * sum(KsiT2)' ),...
             spec.s, spec.s)';

        % Minimize the likelihood function and retrieve the matricex B and
        % Lambdas, BL stores  B first, then it has diagonal entries of Lambdas using
        % function Minimize
        convergenceValues = 0;
        convergenceValues(1,1)=1;
        convergenceValues(1,2)=1;
        innerRound = 1;
        while convergenceValues(innerRound,2) > 0.000000001
            u = GetResiduals(T, y, Z, spec.lags, theta(:, totalRounds));
            [BL, convergenceValues(innerRound+1,1)] = MinimizeN(shoothedRegimes, T, Tm, Lambda, B, spec, u, theta, totalRounds);
            
            convergenceValues(innerRound+1,2) = abs( ( convergenceValues(innerRound+1,1) - convergenceValues(innerRound,1) ) / convergenceValues(innerRound,1) );

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
                T2 = T2 + kron( ( Z' * diag(shoothedRegimes(2:T(1,1)-spec.lags+1,position)) * Z ), GetSigmaForRegime(Sigma, T, position)^-1) ;
                T3 = T3 + kron( Z' * diag(shoothedRegimes(2:T(1,1)-spec.lags+1,position)), GetSigmaForRegime(Sigma, T, position)^-1);
            end
            theta(:, totalRounds+1) = T2^-1 * T3 * reshape(y(spec.lags+1:T(1,1),:)', (T(1,1)-spec.lags) * T(1,2),1) ;
            innerRound = innerRound + 1;
            totalRounds =totalRounds + 1;
        end
        % update initial p
        filteredProbabilities(1,:) = shoothedRegimes(1,:);   

        fprintf(1,'Current round of iteration of EM is ');
        fprintf(1,'%d \n', outerRound);

        outerRound = outerRound + 1;
        
    end
    
thetaFinal = theta(:, totalRounds);
likelihoodFinal = [logValue(outerRound,2), convergenceValues(innerRound,1)];