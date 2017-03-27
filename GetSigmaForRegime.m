% Retrieve the covariance matrix for prespecified position
function[Sigma] = GetSigmaForRegime(input_sigma, T, position)

Sigma = input_sigma (1 + T(1,2)*(position-1) : T(1,2) * position, :) ;

