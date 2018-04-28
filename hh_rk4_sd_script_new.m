function x1 = hh_rk4_sd_script_new (v,I,c)
% v = -60; 
% an = 0.01*(v+50)/(1-exp(-(v+50)/10));
% bn = 0.125*exp(-(v+60)/80);
% am = 0.1*(v+35)/(1-exp(-(v+35)/10));
% bm = 4*exp(-(v+60)/18);
% ah = 0.07*exp(-(v+60)/20);
% bh = 1/(1+exp(-(v+30)/10));

an = c(7)./(1+exp(c(8).*(v - c(9))));
bn = c(10)./(1+exp(c(11).*(v - c(12))));
am = c(13)./(1+exp(c(14).*(v - c(15))));
bm = c(16)./(1+exp(c(17).*(v - c(18))));
ah = c(19)./(1+exp(c(20).*(v - c(21))));
bh = c(22)./(1+exp(c(23).*(v - c(24))));

n = an/(an+bn); m = am/(am+bm); h = ah/(ah+bh);
t = 0:0.0126:126; 
y(1,:) = [n m h v];
for i=1:(length(t)-2)    
    ytmp1 = hh_rk4_new(I(i),t(i+2)-t(i), y(i,:)',c);
    ytmp2 = hh_rk4_new(I(i),t(i+1)-t(i), y(i,:)',c);
    ytmp3 = hh_rk4_new(I(i),t(i+1)-t(i), ytmp2,c);    
    if ((16/15)*(ytmp3 - ytmp1) <= (5e-3)*(ytmp3))
        y(i+1,:) = ytmp1';
    else
        y(i+1,:) = ytmp2';        
    end
end
ytmp = hh_rk4_new(I(i),t(length(t))-t(length(t)-1), y(length(t)-1,:)',c);
y(length(t),:) = ytmp';
figure;
plot(t,y(:,1)); hold on;
plot(t,y(:,2)); hold on;
plot(t,y(:,3)); grid on;
figure;
plot(t,y(:,4)); grid on
x1 = y(:,4);
x1d = [diff(x1);0];
figure; plot(x1,x1d);grid on
end