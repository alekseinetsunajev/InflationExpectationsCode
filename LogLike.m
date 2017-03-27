function logLikelihood = LogLike(x, KsiT,  Tm, spec, u, T )
logLikelihood = 0;
TMP = T(1,2);
for position = 1 : spec.s
    A(1 + T(1,2)*(position-1) : T(1,2) * position, :) = (repmat(KsiT(2:T(1,1) - spec.lags + 1, position),1, T(1,2) ) .* u )' * u;
end

logLikelihood = logLikelihood + (T(1,1)-spec.lags)*log( abs( det( x(:, 1: T(1,2)) ) ) ) + ...
0.5*trace( x(:, 1: T(1,2))'^-1 * x(:, 1:T(1,2))^-1 * GetSigmaForRegime(A, T, 1));

for position = 2:spec.s
    logLikelihood = logLikelihood + Tm(1,position)/2 * log( det( diag( x(:,TMP+position-1) ) ) ) + ...
    0.5*trace( x(:,1: TMP)'^-1 * diag( x(:, TMP+position-1) )^-1 * x(:, 1: TMP)^-1 * GetSigmaForRegime(A, T, position) );
end