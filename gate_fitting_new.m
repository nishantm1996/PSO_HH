function [] = gate_fitting_new(a,b)
% v = -100:0.001:100;
% an = 0.01.*(v+50)./(1-exp(-(v+50)./10));
% bn = 0.125.*exp(-(v+60)./80);
% am = 0.1.*(v+35)./(1-exp(-(v+35)./10));
% bm = 4.*exp(-(v+60)./18);
% ah = 0.07.*exp(-(v+60)./20);
% bh = 1./(1+exp(-(v+30)./10));

an_ = [1.5 -0.03 25];
bn_ = [0.22 0.018 -52];

am_ = [15.75 -0.03 45];
bm_ = [40 0.075 -90];

ah_ = [0.57 0.065 -90];
bh_ = [1 -0.1 -30];

% figure;
v = min(a):max(a);
anew = an_(1)./(1+exp(an_(2).*(v - an_(3))));
subplot(3,2,1);
plot(anew); hold on
title('an');
v = min(b):max(b);
anew = an_(1)./(1+exp(an_(2).*(v - an_(3))));
plot(anew); grid on
% plot(an); hold on; 
% legend('Original','PSO');
% 
% figure;
v = min(a):max(a);
bnew = bn_(1)./(1+exp(bn_(2).*(v - bn_(3))));
subplot(3,2,2);
plot(bnew); hold on
title('bn');
v = min(b):max(b);
bnew = bn_(1)./(1+exp(bn_(2).*(v - bn_(3))));
plot(bnew); grid on
% plot(bn); hold on; 
% legend('Original','PSO');
% 
% figure;
v = min(a):max(a);
amew = am_(1)./(1+exp(am_(2).*(v - am_(3))));
subplot(3,2,3);
plot(amew); hold on
title('am');
v = min(b):max(b);
amew = am_(1)./(1+exp(am_(2).*(v - am_(3))));
plot(amew); grid on
% plot(am); hold on; 
% legend('Original','PSO');
% 
% figure;
v = min(a):max(a);
bmew = bm_(1)./(1+exp(bm_(2).*(v - bm_(3))));
subplot(3,2,4);
plot(bmew); hold on
title('bm');
v = min(b):max(b);
bmew = bm_(1)./(1+exp(bm_(2).*(v - bm_(3))));
plot(bmew); grid on
% plot(bm); hold on; 
% legend('Original','PSO');
% 
% figure;
v = min(a):max(a);
ahew = ah_(1)./(1+exp(ah_(2).*(v - ah_(3))));
subplot(3,2,5);
plot(ahew); hold on
title('ah');
v = min(b):max(b);
ahew = ah_(1)./(1+exp(ah_(2).*(v - ah_(3))));
plot(ahew); grid on
% plot(ah); hold on; 
% legend('Original','PSO');
% 
% figure;
v = min(a):max(a);
bhew = bh_(1)./(1+exp(bh_(2).*(v - bh_(3))));
subplot(3,2,6);
plot(bhew); hold on
title('bh');
v = min(b):max(b);
bhew = bh_(1)./(1+exp(bh_(2).*(v - bh_(3))));
plot(bhew); grid on
% plot(bh); hold on; 
% legend('Original','PSO');
end


