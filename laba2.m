clear; clc;

A = 3; % начальные условия
I0 = A-1;
D0 = 1;
tp = 60;
fa = 0.8834; % 10^9
ts = 10/fa;
te = 1000;
num = 50;
famin = fa*0.05;
famax = fa*1.25;
step = (famax-famin)/num;

c=1;
u = zeros(num+1, 3);
for fa_curr = famin:step:famax
  [~, yy] = auxiliary(A, I0, D0, tp, ts, te, fa_curr);
  u(c, 1) = fa_curr;
  u(c, 2) = min(yy(:,1));
  u(c, 3) = max(yy(:,1));
  c=c+1;

end
figure;
plot(u(:,1), u(:,2), 'b-o', u(:,1), u(:,3), 'r-o');
legend('Imin','Imax');
xlabel('f, 10^9')
ylabel('I(f)')