function res=cycloidal(t,T,S0,dS)
% cycloidal motion curve
res.pos=S0+dS*(t/T-sin(2*pi*t/T)/(2*pi));
res.vel=dS*(1-cos(2*pi*t/T))/T;
res.acc=2*pi*dS/T^2*sin(2*pi*t/T);
end