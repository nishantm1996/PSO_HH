function yn1 = hh_rk4_new (I,h,yn,c)
k1 = coeff_a_b_rk4_new(0,I,yn,c);
k2 = coeff_a_b_rk4_new(0,I,yn + k1.*(h/2),c);
k3 = coeff_a_b_rk4_new(0,I,yn + k2.*(h/2),c);
k4 = coeff_a_b_rk4_new(0,I,yn + k3.*h,c);
yn1 = yn + h.*(k1 + 2*k2 + 2*k3 + k4)./6;
end