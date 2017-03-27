function[C] = GetCoefficient(Theta, T, Position)

cnt1 = (Position-1)*T(1,2)^2 + T(1,2) + 1;   
cnt2 = (Position)*T(1,2)^2 + T(1,2);

C = reshape(Theta(cnt1:cnt2), T(1,2), T(1,2) );