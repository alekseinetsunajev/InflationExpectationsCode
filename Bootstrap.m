% This function estimates the Impulse Responses and computes various
% bootstrap Confidence Intervals. Procedure follows Lütkepohl (2005) D.3,
% p.709.
function[] = Bootstrap(spec, VARcoefficients, B, smoothedRegimeProbs, Lambda, T, y, Z, h, brep)

% To compute IRs with actual data:
if B(1,1) < 0
    B(:,1) = -1 * B(:,1);
end
if B(2,2) < 0
    B(:,2) = -1 * B(:,2);
end
[pointEstimateImpulseResponse, pointEstimateAccumulatedImpulseResponse] = CalculateImpulseResponses(spec, VARcoefficients, B,  T, h);

residuals = GetResiduals(T, y, Z, spec.lags, VARcoefficients);
constant = VARcoefficients(1:T(1,2))' ;

%% The bootstrap loop
bootstrapImpulseResponses = zeros(T(1,2) ^ 2,  (h+1), brep);
for n=1:brep
    fprintf(1,'Current bootstrap replication ');
    fprintf(1,'%d \n', n);
     % Set initial values:
    bootstrapData = zeros(T(1,1), T(1,2) );
    bootstrapData(1: spec.lags, :) = y(1:spec.lags, :); % Pre-sample values    
    
    %% Fixed design wild bootstrap series 
    
    % Rademacher distribution:
    rademacherNumbers = 2*round(rand(T(1,1)-spec.lags, 1))-1;
    
    %  use the Rademacher distribution,
    for t=1:(T(1,1)-spec.lags)
        bootstrapResiduals(t, :) = (residuals(t, :)' * rademacherNumbers(t,:) )';
    end
          
    % Generate wild bootstraped series,
    for t = spec.lags+1 : T(1,1)
        if spec.trend == 1
            for i=T(1,2):T(1,2):T(1,2)*spec.lags
                bootstrapData(t, :)  = bootstrapData(t, :) + (GetCoefficient(VARcoefficients, T, i/T(1,2)) * y(t-i/T(1,2), :)')' ;
            end
            bootstrapData(t, :) = bootstrapData(t, :) + constant + t * VARcoefficients(T(1,2)+1:2*T(1,2))'  + bootstrapResiduals(t-spec.lags, :);
        else
            for i=T(1,2):T(1,2):T(1,2)*spec.lags
                bootstrapData(t, :)  = bootstrapData(t, :) + (GetCoefficient(VARcoefficients, T, i/T(1,2)) * y(t-i/T(1,2), :)')' ;
            end
            bootstrapData(t, :) = bootstrapData(t, :) + constant + bootstrapResiduals(t-spec.lags, :);
        end
    end

    % Obtain parameters based on bootstrapped series    
    Y = bootstrapData(spec.lags+1:end, :);
    A1 = 0;
    A2 = 0;
    for i = 1:T(1,1) - spec.lags;
        bootstrappedLaggedValues(i,1) = 1;            % constant
        if spec.trend == 1      % with trend
            bootstrappedLaggedValues(i,2) = i;
            for j=1:spec.lags
                bootstrappedLaggedValues(i,(j-1)*T(1,2)+2+1 : j*T(1,2)+1+1 )= y(spec.lags + i - j , :);  % lagged observations 
            end
        else
             for j=1:spec.lags
                bootstrappedLaggedValues(i,(j-1)*T(1,2)+2 : j*T(1,2)+1 )= y(spec.lags + i - j , :);  % lagged observations 
            end           
        end
    end
    
    roundOfIteration = 1;
    convergenceParameters(1,1)=1;
    convergenceParameters(1,2)=1;
    
    observedRegimes = sum(smoothedRegimeProbs ) - smoothedRegimeProbs(1,:);
    bootstrappedVARcoefficient  = VARcoefficients;
    estimatedResiduals =  GetResiduals(T, bootstrapData, bootstrappedLaggedValues, spec.lags, bootstrappedVARcoefficient);

    for position = 1:spec.s
        Sigma(1 + T(1,2)*(position-1) : T(1,2) * position, : )= B * ...
                            Lambda(1 + T(1,2)*(position-1) : T(1,2) * position, :) * B' ;
    end
    bootstrapVectorisedLambda(:, 1) = reshape(Lambda,  T(1,2)^2*spec.s, 1 );

    while convergenceParameters(roundOfIteration,2) > 0.0000000001
      
       [BL, convergenceParameters(roundOfIteration+1,1)] = MinimizeN(smoothedRegimeProbs, T, observedRegimes, Lambda, B, spec, estimatedResiduals, bootstrappedVARcoefficient, 1);

       % Store B 
       B = BL(:,1: T(1,2));
        
       % store percentage change
       convergenceParameters(roundOfIteration+1,2) = abs( ( convergenceParameters(roundOfIteration+1,1) - convergenceParameters(roundOfIteration,1) ) / convergenceParameters(roundOfIteration,1) );

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
            T2 = T2 + kron( ( bootstrappedLaggedValues' * diag(smoothedRegimeProbs(2:T(1,1)-spec.lags+1,position)) * bootstrappedLaggedValues ), GetSigmaForRegime(Sigma, T, position)^-1) ;
            T3 = T3 + kron( bootstrappedLaggedValues' * diag(smoothedRegimeProbs(2:T(1,1)-spec.lags+1,position)), GetSigmaForRegime(Sigma, T, position)^-1);
       end
       bootstrappedVARcoefficient = T2^-1 * T3 * reshape(Y', (T(1,1)-spec.lags) * T(1,2),1) ;
       estimatedResiduals = GetResiduals(T, bootstrapData, bootstrappedLaggedValues, spec.lags, bootstrappedVARcoefficient);        
       roundOfIteration = roundOfIteration + 1;
    end
    bootsrapVectorisedB(:, n) = reshape(B, T(1,2)^2, 1);
    bootstrapVectorisedLambda(:, n+1) = reshape(Lambda,  T(1,2)^2*spec.s, 1 );    
      
    [bootstrapImpulseResponses(:,:, n), bootstrapAccumulatedImpulseResponses(:,:, n)] = CalculateImpulseResponses(spec, bootstrappedVARcoefficient, B,  T, h);
end

%% Compute Confidence Intervals
q = 0.975;
confIntervalHigh = quantile(bootstrapAccumulatedImpulseResponses,q,3);
confIntervalLow = quantile(bootstrapAccumulatedImpulseResponses,1-q,3);

% rescale the impulse responses s.t. the shocks have a unit impact
responsesToPlot = pointEstimateAccumulatedImpulseResponse;
responsesToPlot([1, 3],:) = pointEstimateImpulseResponse([1, 3],:);
confIntervalHigh([1, 3], :) = quantile(bootstrapImpulseResponses([1, 3],:,:),q,3);
confIntervalLow([1, 3], :) = quantile(bootstrapImpulseResponses([1, 3],:,:),1-q,3);

indexOfPosition = 0;
positon = [1 2 3 4];        % position
for i = 1:T(1,2):T(1,2)^2
    confIntervalLow(i: i + T(1,2) - 1, :) = confIntervalLow(i : i + T(1,2) - 1, :) ./ pointEstimateImpulseResponse(positon(i + indexOfPosition),1);
    confIntervalHigh(i : i + T(1,2) - 1, :) = confIntervalHigh(i : i + T(1,2) - 1, :) ./ pointEstimateImpulseResponse(positon(i + indexOfPosition) ,1);
    responsesToPlot(i:i + T(1,2) - 1, :) = responsesToPlot(i:i + T(1,2) - 1, :) ./ pointEstimateImpulseResponse(positon(i + indexOfPosition),1);
    indexOfPosition = indexOfPosition + 1;
end

%% Plot of Impulse Responses with Confidence Intervals
% PICTURE FOR INFLATION EXPECTATIONS
figure('Name','Impulse Responses','NumberTitle','off')
periods=0:h;
Vec = reshape(reshape((1:T(1,2)^2), T(1,2), T(1,2))', T(1,2)^2, 1)';
for p = 1 : T(1,2)^2
                subplot(T(1,2), T(1,2), p);
                plot(periods, responsesToPlot(Vec(p),:),'k',  'LineWidth', 1.5);
                hold on;
                plot(periods, confIntervalHigh(Vec(p),:),'k--',  'LineWidth', 1.5);
                plot(periods, confIntervalLow(Vec(p),:),'k--' ,  'LineWidth', 1.5);
                plot(periods, zeros(1, h+1), ':' );
                hold off;
                xlim(gca, [0 h]);
end
