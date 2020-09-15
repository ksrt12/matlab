clear; clc;
% part1();
part2();

function part1()
subplot(1,3,1)
fplot(@(x) func(x), [-10 10]);
xlabel('x'); ylabel("y");
subplot(1,3,2)
fplot(@(x) func1(x), [-10 10]);
xlabel('x'); ylabel("y'");
subplot(1,3,3)
fplot(@(x) func2(x), [-10 10]);
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
A = 2; % начальные условия
I0 = 5;
D0 = 1;
tp = 60*10^-12;
tc = 1*10^-9;
ts = 0;
te = 50;
lab_1_task_2(I0, D0, A, tp, tc, ts, te);
end

