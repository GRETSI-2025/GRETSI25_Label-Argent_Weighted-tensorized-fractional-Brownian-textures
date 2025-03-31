% Simulation of Weighted Tensoriel Fractional Brownian Textures (WTFBT)
% 
% Authors: Claire Launay, Béatrice Vedel
% Date: 02/2024 (Modified in 03/2025)
% Associated with the preprint "Weighted tensorized fractional Brownian
% textures", Céline Esser, Claire Launay, Laurent Loosvelt, Béatrice Vedel, March 2025
%

close all
clear

%% Test stationary rectangular increments for 100 images

M = 9;
Mx = 2^M; My = 2^M;
H = 0.3;
alpha = 0.5;
sz = 2^(M-1);

mean_est_incr = zeros(sz,sz);
var_est_incr = zeros(sz,sz);
skew_est_incr = zeros(sz,sz);

max_iter = 100;
mean_rect_incr_all = zeros(1,max_iter);
var_rect_incr_all = zeros(1,max_iter);
skew_rect_incr_all = zeros(1,max_iter);
mean_est = zeros(1,max_iter);
var_est = zeros(1,max_iter);
skew_est = zeros(1,max_iter);

rng(0)

for iter = 1:max_iter
    display(['Iter : ', num2str(iter)])
    z1 = randn(2*My,2*Mx); z2 = randn(2*My,2*Mx);
    z = z1 + 1i*z2;
        
    freqx= -Mx+1:Mx ; freqy = -My+1:My;
    [XX,YY] = meshgrid(pi*freqx,pi*freqy);
    phi = (min(abs(XX),abs(YY)).^((1-alpha)*H+1/2)  .*  max(abs(XX),abs(YY)).^((1+alpha)*H+1/2));
    g = 1./phi;
    g(My,:) = 0; g(:,Mx)=0;
    
    y1 = fft(z.*g, 2*Mx,2);
    y2 = fft(y1 - ones(2*My,1)*y1(My,:),2*My,1);
    y = y2 - y2(:,Mx)*ones(1,2*Mx);
        
    y = pi* y(My:end, Mx:end); 
    x = real(y);

    mean_est(iter) = mean(x(:));
    var_est(iter) = var(x(:));
    skew_est(iter) = skewness(x(:));


    for indi = 1:sz
        for indj = 1:sz
            x_rect_incr = x(indi:sz+indi-1,indj:sz+indj-1) -ones(sz,1)*x(indi,indj:sz+indj-1) - x(indi:sz+indi-1,indj)*ones(1,sz)+x(indi,indj);
    
            mean_est_incr(indi,indj) = mean(x_rect_incr(:));
            var_est_incr(indi,indj) = var(x_rect_incr(:));
            skew_est_incr(indi,indj) = skewness(x_rect_incr(:));
    
        end
    end
    
    mean_rect_incr_all(iter) = mean(mean_est_incr(:));
    var_rect_incr_all(iter) = mean(var_est_incr(:));
    skew_rect_incr_all(iter) = mean(skew_est_incr(:));

end


disp(['Empirical mean is ', num2str(mean(mean_est)), ', empirical variance is ', num2str(mean(var_est)), ', and empirical skewness is ', num2str(mean(skew_est))])
disp(['Rectangular-increments mean is ', num2str(mean(mean_rect_incr_all)), ', rectangular-increments variance is ', num2str(mean(var_rect_incr_all)), ', and rectangular-increments skewness is ', num2str(mean(skew_rect_incr_all))])


%% Test autosimilarity for 100 images

M = 9;
Mx = 2^M; My = 2^M;
H = 0.3;
alpha = 0.5;
a_autosim = [2,3,4,5,6,7,8];
max_iter = 100;

mean_est_autosim = zeros(max_iter,length(a_autosim));
var_est_autosim = zeros(max_iter,length(a_autosim));
skew_est_autosim = zeros(max_iter,length(a_autosim));

max_iter = 100;
mean_autosim_all = zeros(1,max_iter);
var_autosim_all = zeros(1,max_iter);
skew_autosim_all = zeros(1,max_iter);
mean_est = zeros(1,max_iter);
var_est = zeros(1,max_iter);
skew_est = zeros(1,max_iter);

rng(0)

for iter = 1:max_iter
    display(['Iter : ', num2str(iter)])
    z1 = randn(2*My,2*Mx); z2 = randn(2*My,2*Mx);
    z = z1 + 1i*z2;
        
    freqx= -Mx+1:Mx ; freqy = -My+1:My;
    [XX,YY] = meshgrid(pi*freqx,pi*freqy);
    phi = (min(abs(XX),abs(YY)).^((1-alpha)*H+1/2)  .*  max(abs(XX),abs(YY)).^((1+alpha)*H+1/2));
    g = 1./phi;
    g(My,:) = 0; g(:,Mx)=0;
    
    y1 = fft(z.*g, 2*Mx,2);
    y2 = fft(y1 - ones(2*My,1)*y1(My,:),2*My,1);
    y = y2 - y2(:,Mx)*ones(1,2*Mx);
        
    y = pi* y(My+1:end, Mx+1:end);
    x = real(y);

    mean_est(iter) = mean(x(:));
    var_est(iter) = var(x(:));
    skew_est(iter) = skewness(x(:));


    for indi = 1:length(a_autosim)
        x_autosim = x(1:a_autosim(indi):end,1:a_autosim(indi):end);
    
        mean_est_autosim(iter,indi) = (1/a_autosim(indi)^(2*H))*mean(x_autosim(:));
        var_est_autosim(iter,indi) = (1/a_autosim(indi)^(4*H))*var(x_autosim(:));
        skew_est_autosim(iter,indi) = a_autosim(indi)^(6*H) *skewness(x_autosim(:));

    end
    
    mean_autosim_all(iter) = mean(mean_est_autosim(iter,:));
    var_autosim_all(iter) = mean(var_est_autosim(iter,:));
    skew_autosim_all(iter) = mean(skew_est_autosim(iter,:));

end



disp(['Empirical mean is ', num2str(mean(mean_est)), ', empirical variance is ', num2str(mean(var_est)), ', and empirical skewness is ', num2str(mean(skew_est))])
disp(['Rescaled field mean is ', num2str(mean(mean_autosim_all)), ', rescaled field variance is ', num2str(mean(var_autosim_all)), ', and rescaled field skewness is ', num2str(mean(skew_autosim_all))])



