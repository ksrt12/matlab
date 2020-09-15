function lab_1_task_2(I0, D0, A, tp, tc, ts, te)

gam = tp/tc;
opt = odeset('RelTol',0.0001);

[tt,yy]=ode45(@func, [ts te], [I0 D0], opt);

plot(tt*tp/(10^(-9)),yy(:,1),'-o')
xlabel('t (ns)')
ylabel('I(t)')

function dy = func(t, y)
  dI = y(1)*(y(2)-1);

  dD = gam*(A - y(2)*(1+y(1)));

  dy = [dI dD].';
end

end
