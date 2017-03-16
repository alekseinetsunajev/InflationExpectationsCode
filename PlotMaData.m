figure
plot(data(7:end, 1), 'k');
hold on;
plot(yMA(:, 1), 'b:');

figure
plot(data(7:end, 2), 'k');
hold on;
plot(yMA(:, 2), 'b:');

figure
plot( yAccumulated(7:end, 2), 'k');
hold on;
plot(yAccumulatedMA(:, 2), 'b:');
