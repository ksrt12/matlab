clear; clc;
A = 3; % начальные условия

main(A);
main(A*0.75);

function main(A)

D0 = 1;
tp = 60; % пс
fa = 0.054053; % ГГц
ts = 0; % t_start
num = 50; % количество точек
famin = fa*0.05; % 5%
famax = fa*1.25; % 125%
step = (famax-famin)/num;
I0 = A-1;
c = 1;
u = zeros(num+1, 3);

for fa_curr = famin:step:famax
  [~, yy] = auxiliary(A, I0, D0, tp, ts, 10/fa_curr, fa_curr);
  yy((1:110), 1) = I0;
  u(c, 1) = fa_curr;
  u(c, 2) = min(yy(:,1));
  u(c, 3) = max(yy(:,1));
  c = c+1;
end
% writematrix(u,strcat(""+A,'.xlsx'),'Sheet',1); % сохранение в таблицу excel
figure('Name',"A="+A);
plot(u(:,1), u(:,2), 'b-o', u(:,1), u(:,3), 'r-o'); % построение графиков
legend('Imin','Imax');
title("I(f) при A="+A);
xlabel('f, tp^-1')
ylabel('I(f)')
end