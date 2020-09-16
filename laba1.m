clear; clc;
part1();
part2();

function part1()
range = [-15 15];
subplot(1,3,1)
fplot(@(x) func(x), range);
xlabel('x'); ylabel("y");
subplot(1,3,2)
fplot(@(x) func1(x), range);
xlabel('x'); ylabel("y'");
subplot(1,3,3)
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
te = 100;
auxiliary(A, I0, D0, tp, te);
end

