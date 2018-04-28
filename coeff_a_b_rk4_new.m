function d = coeff_a_b_rk4_new (t,I,x,c)
% I = 0.1;
cm = 0.01;
ena = c(1);
ek = c(2);
el = c(3);
gna = c(4);
gk = c(5);
gl = c(6);
d = zeros(4,1);
an = c(7)./(1+exp(c(8).*(x(4) - c(9))));
bn = c(10)./(1+exp(c(11).*(x(4) - c(12))));
am = c(13)./(1+exp(c(14).*(x(4) - c(15))));
bm = c(16)./(1+exp(c(17).*(x(4) - c(18))));
ah = c(19)./(1+exp(c(20).*(x(4) - c(21))));
bh = c(22)./(1+exp(c(23).*(x(4) - c(24))));
d(1) = an*(1-x(1)) - bn*x(1);
d(2) = am*(1-x(2)) - bm*x(2);
d(3) = ah*(1-x(3)) - bh*x(3);
d(4) = (1/cm)*(I - gna*(x(2)^3)*x(3)*(x(4)-ena) - gk*(x(1)^4)*(x(4)-ek) - gl*(x(4)-el));
end