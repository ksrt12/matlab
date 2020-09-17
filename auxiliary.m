function [tt, yy] = auxiliary(A, I0, D0, tp_p, ts, te, fa)

if (nargin > 6)
    sa = A*0.1;
else
    sa = 0;
    fa = pi;
end

range = [ts te];
tp = tp_p*10^-12;
tc = 10^(-9);
gam = tp/tc;
opt = odeset('RelTol',0.0001);

[tt,yy]=ode23(@func, range, [I0 D0], opt);

function dy = func(t, y)
  dI = y(1)*(y(2)-1);

  dD = gam*(A + sa*sin(2*pi*fa*t) - y(2)*(1+y(1)));

  dy = [dI dD].';
end

end

