function lab_1_task_2(A, I0, D0, tp_p, plot_opt)

if nargin < 5
    plot_opt = '';
end

tp = tp_p*10^-12;
tc = 10^(-9);
gam = tp/tc;
opt = odeset('RelTol',0.0001);

ts = 0; te = 50;

[tt,yy]=ode45(@func, [ts te], [I0 D0], opt);

figure;
plot(tt*tp/tc,yy(:,1),plot_opt)
xlabel('t (ns)')
ylabel('I(t)')

function dy = func(~, y)
  dI = y(1)*(y(2)-1);

  dD = gam*(A - y(2)*(1+y(1)));

  dy = [dI dD].';
end

end
