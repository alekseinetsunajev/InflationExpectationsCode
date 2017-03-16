function residuals = GetResiduals(T, y, Z, lags, Theta)

residuals = y(lags + 1 : T(1,1), :) - Z * reshape(Theta, T(1,2),  lags*T(1,2) + 1)'; 