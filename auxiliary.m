function auxiliary(A, I0, D0, tp_p, te, laba)

tp = tp_p*10^-12;
tc = 10^(-9);
gam = tp/tc;
opt = odeset('RelTol',0.0001);

ts = 0;

[tt,yy]=ode45(@func, [ts te], [I0 D0], opt);

figure;
plot(tt*tp/tc,yy(:,1),'-o')
xlabel('t (ns)')
ylabel('I(t)')

function dy = func(~, y)
  dI = y(1)*(y(2)-1);

  dD = gam*(A - y(2)*(1+y(1)));

  dy = [dI dD].';
end

end
