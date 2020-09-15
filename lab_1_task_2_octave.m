I0 = 2;
D0 = 2;

A = 3;

tau_p = 200*10^(-12);

tau_c = 1*10^(-9);

gam = tau_p/tau_c;

function dy = model(t,y)
 
  dI = y(1)*(y(2)-1);  

  dD = gam*(A - y(2)*(1+y(1)));
  
  dy = [dI dD];
  
endfunction

tstart = 0;
tend = 50;

opt=odeset('RelTol',0.0001);

[tt,yy]=ode23(@model, [tstart, tend], [I0; D0], opt);

plot(tt*tau_p/(10^(-9)),yy(:,1))
xlabel('t (ns)')
ylabel('I(t)')
end