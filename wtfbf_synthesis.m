% Simulation of Weighted Tensoriel Fractional Brownian Textures (WTFBT)
% 
% Authors: Claire Launay, Béatrice Vedel
% Date: 02/2024 (Modified in 03/2025)
% Associated with the preprint "Weighted tensorized fractional Brownian
% textures", Céline Esser, Claire Launay, Laurent Loosveldt, Béatrice Vedel, March 2025
%

close all
clear

%% General case of WTFBF of size (Mx+1)x(My+1), with parameters H and alpha 

M = 10;
Mx = 2^M; My = 2^M;

for H = [0.3,0.7]
    for alpha = [0,0.5,1]
        rng(3)
        
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
        
        figure; 
        imagesc(x(1:2:end,1:2:end)); 
        colormap('gray');
        title(['WTFBT for H = ', num2str(H), ' and alpha = ', num2str(alpha)]); 
        axis equal; axis off
    
    end
end

%% Fractional Brownian sheet (alpha = 0)

H = 0.3;
M = 10;
Mx = 2^M; My = 2^M;

rng(3)

z1 = randn(2*My,2*Mx); z2 = randn(2*My,2*Mx);
z = z1 + 1i*z2;

freqx= -Mx+1:Mx ; freqy = -My+1:My;
[XX,YY] = meshgrid(pi*freqx,pi*freqy);
phi = (abs(XX).^(H+1/2).* abs(YY).^(H+1/2));
g = 1./phi;
g(My,:) = 0; g(:,Mx)=0;

y1 = fft(z.*g, 2*Mx,2);
y2 = fft(y1 - ones(2*My,1)*y1(My,:),2*My,1);
y = y2 - y2(:,Mx)*ones(1,2*Mx);

y = pi* y(My:end, Mx:end);

x = real(y);

figure; 
imagesc(x(1:2:end,1:2:end)); 
colormap('gray');
title(['Fractional Brownian sheet for H = ', num2str(H),]); 
axis equal; axis off



%% Fractional Brownian field (close to alpha = 1)

H = 0.3;
M = 10;
Mx = 2^M; My = 2^M;

rng(3)

z1 = randn(2*My,2*Mx); z2 = randn(2*My,2*Mx);
z = z1 + 1i*z2;

freqx= -Mx+1:Mx ; freqy = -My+1:My;
[XX,YY] = meshgrid(pi*freqx,pi*freqy);
phi = sqrt(abs(XX).^2 + abs(YY).^2).^(H+1);
g = 1./phi;
g(My,Mx)=0;

y = pi* fft2(z.*g);
y = y(My:end, Mx:end);
x = real(y-y(1,1));

figure; 
imagesc(x(1:2:end,1:2:end)); 
colormap('gray'); 
title(['Fractional Brownian field for H = ', num2str(H)]); axis equal; axis off



%% Anisotropic WTFBT with beta1 != beta2 
%close all;

H_all = [0.4,0.6];
beta1_all = [0.7,0.85];
beta2_all = [1.3,1.15];
M = 10;
Mx = 2^M; My = 2^M;


for k=1:2
    H = H_all(k);
    beta1 = beta1_all(k);
    beta2 = beta2_all(k);

    if (max(beta1,beta2)-1>2*H)||(2*H >3*min(beta1,beta2)-1)
        error('Constraints on the parameters ( max(beta1,beta2)-1 < 2H < 3min(beta1,beta2)-1 ) are not respected.')
    end
    
    for alpha = [0,0.5,1]
        rng(3)
        
        z1 = randn(2*My,2*Mx); z2 = randn(2*My,2*Mx);
        z = z1 + 1i*z2;
        
        freqx= -Mx+1:Mx; freqy = -My+1:My;
        [XX,YY] = meshgrid(pi*freqx,pi*freqy);
        phi = (min(abs(XX).^(1/beta1),abs(YY).^(1/beta2)).^((1-alpha)*H+1/2)  .*  max(abs(XX).^(1/beta1),abs(YY).^(1/beta2)).^((1+alpha)*H+1/2));
        g = 1./phi;
        g(My,:) = 0; g(:,Mx)=0;
        y1 = fft(z.*g, 2*Mx,2);
        y2 = fft(y1 - ones(2*My,1)*y1(My,:),2*My,1);
        y = y2 - y2(:,Mx)*ones(1,2*Mx);
        
        y = pi* y(My:end, Mx:end);
        x = real(y);
        
        figure; 
        imagesc(x(1:2:end,1:2:end)); 
        colormap('gray');
        title(['WTFBT for H = ', num2str(H), ', beta1 = ', num2str(beta1), ', beta2 = ', num2str(beta2),  ' and alpha = ', num2str(alpha)]);  axis equal; axis off
    end
end

