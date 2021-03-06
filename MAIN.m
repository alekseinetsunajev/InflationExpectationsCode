%% This is the main program of the code to replicate results of the paper "The Anchoring of Ination Expectations in the Short and in the Long Run "
% Authors: 
% Dieter Nautz
% Aleksei Netsunajev
% Till Strohsal

% Code written by Aleksei Netsunajev, aleksei.netsunajev@fu-berlin.de
clear
fprintf(1,'\n');
fprintf(1,'This is the code to replicate result of the paper "The Anchoring of Inflation Expectations in the Short and in the Long Run "\n');
fprintf(1,'\n');

%% Read in some data (levels-difference bivariate specification)
y = xlsread('data_l_d.xls');
% T stores time horizon and number of variables of the system
T = size(y);

%% Specify the model to be estimated
% spec is a structure, that includes the variables describing the system
spec.s = 2;                 % states
spec.lags = 6;            % lags
spec.trend =  0;         % trend
spec.BQrestrict = 0 ;
% starting value
LM = [0.5];
spec.iterations = numel(LM);
spec.const = ones(spec.s,1);            % s dimentional vector of ones

fprintf(1,'The system has ');
fprintf(1,'%d ', spec.s);
fprintf(1,'states ');
fprintf(1,'%d ', T(1,2));
fprintf(1,'variables and ');
fprintf(1,'%d ', spec.lags);
fprintf(1,'lags\n');

%% DEFINITION OF VARIABLES FOR THE EM ALGORITHM

% Estimation is a structure that holds the estmated parameters that correspond to maximum of log-likelihood:
%P                     Markov matrix
%Theta              Parameter vector
%B                     B matrix of covariance matrix decomposition
%Lambda         Lambda matrix of covariance matrix decomposition
%logL               Log likelihood
%Sigma            covariance matrices, stacked each under another
%KsiT               smoothed probabilities with timing 0,1,2...T-1
%FEVD           forecast error variance decompositions

%% Calculate the matrix Z of a constant and lagged observations
tolValue = 0.0000000001;                              % tolerance value for EM convergence
[Z, A1, A2] = CalculateMatrixOfLaggedValues (y, spec.lags, 0);

%% Estimate models, estimate() is the estimation function, results are
% stored in the Estimation structure
% Form matrix of random numbers
for iteration = 1:spec.iterations
    randNumbers(:,1+ T(1,2)*(iteration -1) :T(1,2)*iteration) = 1/100000*randn(T(1,2),T(1,2));
end

%% Unrestricted model
reply = 0;
fprintf(1,'Estimating model with state invariant B... \n');
 [EstimationUR.Theta(:, iteration), EstimationUR.Sigma(:, 1+ T(1,2)*(iteration -1) :T(1,2)*iteration  ), EstimationUR.KsiT( : , 1+ spec.s*(iteration -1) : spec.s*iteration ), ...
     EstimationUR.P(:, 1 + spec.s*(iteration - 1) : spec.s*iteration ),  EstimationUR.IterationLik(iteration, 1:2), EstimationUR.B(:,1 + T(1,2)*(iteration - 1) : T(1,2)*iteration), ...
     EstimationUR.Lambda(:, 1 + T(1,2)*(iteration - 1) : T(1,2)*iteration), EstimationUR.Eta( : , 1+ spec.s*(iteration -1) : spec.s*iteration ), ...
     EstimationUR.Ksi( : , 1+ spec.s*(iteration -1) : spec.s*iteration )] = ...
            estimate((A1)^-1 * A2, T, y, Z, spec, tolValue, randNumbers(:,1+ T(1,2)*(iteration -1) :T(1,2)*iteration), LM(iteration) );
fprintf(1,'Calculating standard errors... \n');
EstimationUR.logL = EstimationUR.IterationLik(1,1);
[EstimationUR.StErr, EstimationUR.StErr_all] = CalculateStrdErrors(T, EstimationUR, spec, y, Z);
EstimationUR.B_stdErr = results(EstimationUR, spec, T, reply); 

EstimationUR.FEVD1= FEVD(EstimationUR.Theta, EstimationUR.B,  EstimationUR.Lambda, spec, T, 100, 1, 0);
EstimationUR.FEVD2 = FEVD(EstimationUR.Theta, EstimationUR.B,  EstimationUR.Lambda, spec, T, 100, 2, 0);

%% Restricted model
spec.BQrestrict = 1 ;
LM = [0.5];
fprintf(1,'Estimating restricted model... \n');
 [Estimation.Theta(:, iteration), Estimation.Sigma(:, 1+ T(1,2)*(iteration -1) :T(1,2)*iteration  ), Estimation.KsiT( : , 1+ spec.s*(iteration -1) : spec.s*iteration ), ...
     Estimation.P(:, 1 + spec.s*(iteration - 1) : spec.s*iteration ),  Estimation.IterationLik(iteration, 1:2), Estimation.B(:,1 + T(1,2)*(iteration - 1) : T(1,2)*iteration), ...
     Estimation.Lambda(:, 1 + T(1,2)*(iteration - 1) : T(1,2)*iteration), Estimation.Eta( : , 1+ spec.s*(iteration -1) : spec.s*iteration ), ...
     Estimation.Ksi( : , 1+ spec.s*(iteration -1) : spec.s*iteration )] = ...
            estimate((A1)^-1 * A2, T, y, Z, spec, tolValue, randNumbers(:,1+ T(1,2)*(iteration -1) :T(1,2)*iteration), LM(iteration) );

fprintf(1,'Calculating standard errors... \n');
Estimation.logL = Estimation.IterationLik(1,1);
[Estimation.StErr, Estimation.StErr_all] = CalculateStrdErrors(T, Estimation, spec, y, Z);
Estimation.B_stdErr = results(Estimation, spec, T, reply); 

%% calculate FEVD
 Estimation.FEVD1= FEVD(Estimation.Theta, Estimation.B,  Estimation.Lambda, spec, T, 500, 1, 0);
 Estimation.FEVD2 = FEVD(Estimation.Theta, Estimation.B,  Estimation.Lambda, spec, T, 100, 2, 0);

%% Produce impulse responces, last input argument is number of replication
fprintf(1,'Bootstrapping...');
Bootstrap(spec, Estimation.Theta, Estimation.B, Estimation.KsiT, Estimation.Lambda, T, y, Z, 36, 10);

% plot some figures
PlotData();
HistoricalDecomposition(Estimation, T, y, Z, spec);
CalculateAutocorrelation(T, y, Z, spec, Estimation, 24);