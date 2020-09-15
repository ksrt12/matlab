I0 = 2;
D0 = 2;

global A;
A = 3;

tau_p = 200*10^(-12);

tau_c = 1*10^(-9);

global gam;
gam = tau_p/tau_c;

tstart = 0;
tend = 50;

opt=odeset('RelTol',0.0001);

[tt,yy]=ode23(@func, [tstart tend], [I0 D0], opt);

plot(tt*tau_p/(10^(-9)),yy(:,1))
xlabel('t (ns)')
ylabel('I(t)')

function dy = func(t, y)
  global A;
  global gam;
  
  dI = y(1)*(y(2)-1);

  dD = gam*(A - y(2)*(1+y(1)));

  dy = [dI dD].';
end