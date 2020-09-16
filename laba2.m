clear; clc;

A = 3; % начальные условия
I0 = A-1;
D0 = 1;
tp = 60;
fa = 883.4*10^6;
% te = 10/fa;
te=200;
[x,y] = auxiliary(A, I0, D0, tp, te, fa);

figure;
plot(x*tp*10^-3,y(:,1))
xlabel('t (ns)')
ylabel('I(t)')

