function [optim,w,mmg] = pso_script_newer (I) 
tic;
n = 20; c1 = 1.5; c2 = 1.5; w = 0.5+0.5*rand(1,100); itrmax = 100; gmin = 0.01; 
% a = [100 -100 -100 2.5 0.75 0.01 3 -0.1 50 0.5 0.05 -100 25 -0.1 100 80 0.02 -200 1 0.02 -200 2 -0.2 -100];
ptheo = [55.17,-72.14,-49.42,1.2,0.36,0.003,1.5,-0.03,25,0.22,0.018,-52,15.75,-0.03,45,40,0.075,-90,0.57,0.065,-90,1,-0.1,-30];
vtheo = hh_rk4_script_new(-60,I,ptheo);

% x(1,:,:) = [100*rand(1,1,n) -100*rand(1,1,n) -100*rand(1,1,n) 2.5*rand(1,1,n) 0.75*rand(1,1,n) 0.01*rand(1,1,n) 3*rand(1,1,n) -0.1*rand(1,1,n) 50*rand(1,1,n) 0.5*rand(1,1,n) 0.05*rand(1,1,n) -100*rand(1,1,n) 25*rand(1,1,n) -0.1*rand(1,1,n) 100*rand(1,1,n) 80*rand(1,1,n) 0.02*rand(1,1,n) -200*rand(1,1,n) rand(1,1,n) 0.02*rand(1,1,n) -200*rand(1,1,n) 2*rand(1,1,n) -0.2*rand(1,1,n) -100*rand(1,1,n)];
% u(1,:,:) = 0.1.*[100*rand(1,1,n) -100*rand(1,1,n) -100*rand(1,1,n) 2.5*rand(1,1,n) 0.75*rand(1,1,n) 0.01*rand(1,1,n) 3*rand(1,1,n) -0.1*rand(1,1,n) 50*rand(1,1,n) 0.5*rand(1,1,n) 0.05*rand(1,1,n) -100*rand(1,1,n) 25*rand(1,1,n) -0.1*rand(1,1,n) 100*rand(1,1,n) 80*rand(1,1,n) 0.02*rand(1,1,n) -200*rand(1,1,n) rand(1,1,n) 0.02*rand(1,1,n) -200*rand(1,1,n) 2*rand(1,1,n) -0.2*rand(1,1,n) -100*rand(1,1,n)];

for i=1:n
    x(1,:,i) = 0.85.*ptheo + 0.3.*ptheo.*rand(1,24);
    u(1,:,i) = 0.085.*ptheo + 0.03.*ptheo.*rand(1,24);
end

P(1,:,:) = x(1,:,:); optim = zeros(1,24);
g = ones(itrmax,n); gG = ones(1,itrmax);
for i=1:n
    tmp = hh_rk4_script_new(-60,I,x(1,:,i));
    tmp(isnan(tmp)) = 0;
%     g(1,i) = sqrt(mean([(tmp(1:80) - vtheo(1:80)).*(tmp(1:80) - vtheo(1:80));2.*(tmp(81:139) - vtheo(81:139)).*(tmp(81:139) - vtheo(81:139));(tmp(140:1001) - vtheo(140:1001)).*(tmp(140:1001) - vtheo(140:1001))]))/sqrt(mean([vtheo(1:80).*vtheo(1:80);2.*vtheo(81:139).*vtheo(81:139);vtheo(140:1001).*vtheo(140:1001)]));
    g(1,i) = rms(tmp-vtheo)/rms(vtheo);
    if (isnan(g(1,i)))
       g(1,i) = 1;
    end
end
[~,ind] = min(g(1,:));
G(1,:) = P(1,:,ind);

for j=2:itrmax
    for i=1:n
        u(j,:,i) = w(j).*u(j-1,:,i) + (c1*rand).*(P(j-1,:,i) - x(j-1,:,i)) + (c2*rand).*(G(j-1,:) - x(j-1,:,i));
        x(j,:,i) = x(j-1,:,i) + u(j,:,i);
        
%         tmpx = x(j,:,i).*a - x(j,:,i).*x(j,:,i);
%         if (sum(logical(tmpx < 0)) > 0)
%             u(j,logical(tmpx < 0),i) = (0.1.*a(logical(tmpx < 0))).*rand(1,sum(logical(tmpx < 0)));
%             x(j,logical(tmpx < 0),i) = a(logical(tmpx < 0)).*rand(1,sum(logical(tmpx < 0)));            
%         end
        
        tmp = hh_rk4_script_new(-60,I,x(j,:,i));
        tmp(isnan(tmp)) = 0;
%         g(j,i) = sqrt(mean([(tmp(1:80) - vtheo(1:80)).*(tmp(1:80) - vtheo(1:80));2.*(tmp(81:139) - vtheo(81:139)).*(tmp(81:139) - vtheo(81:139));(tmp(140:1001) - vtheo(140:1001)).*(tmp(140:1001) - vtheo(140:1001))]))/rms(vtheo);
%         g(j,i) = sqrt(mean([(tmp(1:80) - vtheo(1:80)).*(tmp(1:80) - vtheo(1:80));2.*(tmp(81:139) - vtheo(81:139)).*(tmp(81:139) - vtheo(81:139));(tmp(140:1001) - vtheo(140:1001)).*(tmp(140:1001) - vtheo(140:1001))]))/sqrt(mean([vtheo(1:80).*vtheo(1:80);2.*vtheo(81:139).*vtheo(81:139);vtheo(140:1001).*vtheo(140:1001)]));
        g(j,i) = rms(tmp-vtheo)/rms(vtheo);
        if (isnan(g(j,i)))
            g(j,i) = 1;
        end
        if (g(j,i) < g(j-1,i))
            P(j,:,i) = x(j,:,i);
        else
            P(j,:,i) = P(j-1,:,i);
            g(j,i) = g(j-1,i);
        end
    end
    [~,ind] = min(g(j,:));
    G(j,:) = P(j,:,ind);
    tmpG = hh_rk4_script_new(-60,I,G(j,:));
    tmpG(isnan(tmpG)) = 0;
%     gG(j) = sqrt(mean([(tmpG(1:80) - vtheo(1:80)).*(tmpG(1:80) - vtheo(1:80));2.*(tmpG(81:139) - vtheo(81:139)).*(tmpG(81:139) - vtheo(81:139));(tmpG(140:1001) - vtheo(140:1001)).*(tmpG(140:1001) - vtheo(140:1001))]))/rms(vtheo);
%     gG(j) = sqrt(mean([(tmpG(1:80) - vtheo(1:80)).*(tmpG(1:80) - vtheo(1:80));2.*(tmpG(81:139) - vtheo(81:139)).*(tmpG(81:139) - vtheo(81:139));(tmpG(140:1001) - vtheo(140:1001)).*(tmpG(140:1001) - vtheo(140:1001))]))/sqrt(mean([vtheo(1:80).*vtheo(1:80);2.*vtheo(81:139).*vtheo(81:139);vtheo(140:1001).*vtheo(140:1001)]));
    gG(j) = rms(tmpG-vtheo)/rms(vtheo);
    if (j>3)
        w(j+1) = exp(-mean(abs(g(j-1,:)-gG(j-1)))./mean(abs(g(j-2,:)-gG(j-2))));
        if (isnan(w(j+1)))
            w(j+1) = 0;
        end
    end
    if (min(g(j,:)) < gmin)
        optim = G(j,:);
        break;
    end
end

if (sum(optim == zeros(1,24)) == 24)
    optim = G(itrmax,:);
end
mmg = min(min(g));
toc;
end

        
            
