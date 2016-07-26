% To obtain the Impulse Responses

function [A, Theta_SR, J, Theta_longrun]=irs(spec, Theta, B,  T, h)
% rescale B such that it has ones on the main diagonal
%B = B * diag(ones(T(1,2), 1) ./ diag(B));

A = zeros(T(1, 2)* (spec.lags) ); % Define the A matrix

for i = 1:spec.lags
      A(1:T(1,2), T(1,2)*(i-1)+1 : T(1,2)*i ) = get_coefficient(Theta, T, i, spec);
end

for m=1:spec.lags-1
    A(m*T(1,2)+1 : m*T(1,2)+T(1,2), T(1,2)*m-T(1,2)+1 : T(1,2) * m) = eye(T(1,2));
end

J = [eye(T(1,2)) zeros(T(1,2), T(1,2)*spec.lags - T(1,2) )  ]; % Define the J matrix

Theta_SR = reshape(J*A^0*J'*B, T(1,2)^2, 1); % Impulse response matrix

for i=1:h
	Theta_SR=([Theta_SR reshape(J*A^i*J'*B, T(1,2)^2, 1)]);
end;

Theta_longrun = Theta_SR(:,1);

for i = 1 : h
    Theta_longrun(:,i+1) =  Theta_longrun(:, i) + Theta_SR(:, i+1);
end
