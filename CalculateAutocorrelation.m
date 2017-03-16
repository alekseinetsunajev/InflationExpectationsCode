function[] = CalculateAutocorrelation(T, y, Z, spec, Estimation, lags)

residuals = GetResiduals(T, y, Z, spec.lags, Estimation.Theta);
for currentVariable = 1: T(1,2)
    [ACF(:, currentVariable), ~, bounds] = autocorr(residuals(:, currentVariable), lags);
end

figure;
subplot(2, 1, 1)
autocorr(residuals(:, 1), lags);
set(gca, 'XLim', [1 lags], 'fontsize', 12);
legend('Sample autocorrelation function','95% confidence intervals')
subplot(2, 1, 2)
autocorr(residuals(:, 2), lags);
set(gca, 'XLim', [1 lags], 'fontsize', 12);

for i = 1:1:T(1,2)
    [~, p, Q, ~] = lbqtest(residuals(:, i), 'lags', lags);
    disp(p);
    disp(Q);
end
