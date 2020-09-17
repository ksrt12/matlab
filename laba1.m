clear; clc;
% part1();
part2();

function part1()
range = [-15 15];
% subplot(1,3,1)
figure;
fplot(@(x) func(x), range);
xlabel('x'); ylabel("y");
% subplot(1,3,2)
figure;
fplot(@(x) func1(x), [-50 50]);
xlabel('x'); ylabel("y'");
% subplot(1,3,3)
figure;
fplot(@(x) func2(x), range);
xlabel('x'); ylabel("y''");
end

function [f] = func(x)
f = x/(sqrt(1+x.^2)); % функция
end

function [f1] = func1(x)
f1 = (-x/((x.^2+1).^1.5)); % 1я производная
end

function [f1] = func2(x)
f1 = (3/((x.^2+1).^2.5)+(x.^2+1).^-1.5); % 2я производная
end

function part2()
A = 3; % начальные условия
I0 = 5;
D0 = 1;
tp = 60;
ts = 0;
te = 100;
[x,y] = auxiliary(A, I0, D0, tp, ts, te);
figure;
plot(x*tp*10^-3,y(:,1),'-o')
xlabel('t (ns)')
ylabel('I(t)')
end

