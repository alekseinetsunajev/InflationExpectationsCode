% To obtain the Impulse Responses
function [responsesShortRun, responsesLongRun] = CalculateImpulseResponses(spec, Theta, B,  T, h)

A = zeros(T(1, 2)* (spec.lags) ); % Define the A matrix

for i = 1:spec.lags
      A(1:T(1,2), T(1,2)*(i-1)+1 : T(1,2)*i ) = GetCoefficient(Theta, T, i);
end

for m=1:spec.lags-1
    A(m*T(1,2)+1 : m*T(1,2)+T(1,2), T(1,2)*m-T(1,2)+1 : T(1,2) * m) = eye(T(1,2));
end

J = [eye(T(1,2)) zeros(T(1,2), T(1,2)*spec.lags - T(1,2) )  ]; % Define the J matrix

responsesShortRun = reshape(J*A^0*J'*B, T(1,2)^2, 1); % Impulse response matrix

for i=1:h
	responsesShortRun=([responsesShortRun reshape(J*A^i*J'*B, T(1,2)^2, 1)]);
end;

responsesLongRun = responsesShortRun(:,1);

for i = 1 : h
    responsesLongRun(:,i+1) =  responsesLongRun(:, i) + responsesShortRun(:, i+1);
end
