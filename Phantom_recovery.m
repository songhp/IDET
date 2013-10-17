%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Compressive Sensing Image Reconstruction                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Coded by Kun Qiu
%Last updated at Feb. 17, 2010

function [PSNR, Count, t] = Phantom_recovery(imgsize, thresh, RadialNum)
% thresh=1e-4;           %Convergence tolerance
%  RadialNum=44; % 38:48 %radial lines 
%
%

% clear all
% clc
% close all

mytime = datestr(now);
disp(sprintf('start at [%s]',mytime));


path(path, './subfunctions');
% path(path,'./Algorithms');

% %Reconstruction
% thresh=1e-4;           %Convergence tolerance for Hard thresholding methods
% 
% RadialNum=44; % 38:48

% imgsize = 128;
Img2D=phantom(imgsize);
img_name='Phantom';
[my,mx]=size(Img2D);
m=my*mx;
Img2Dbar=mean(Img2D(:));
Img2Dcntr=Img2D-Img2Dbar;
Imgcntr=Img2Dcntr(:);
scale=max(Imgcntr)-min(Imgcntr);


%Sampling operator
% number of radial lines in the Fourier domain
[M,Mh,mh,mhi]=LineMask(RadialNum,mx);
OMEGA=mhi;
smpl_pattern=fftshift(M);

Phi=@(z) A_fhp(z,OMEGA);
Phit=@(z) At_fhp(z,OMEGA,mx);

% taking measurements
y=Phi(Imgcntr);
N=length(y);
Nratio=N/m;
%Effective sensing operator
dwt_L=6;                                     %levels of wavelet transform
wav=daubcqf(2);
H=@(z) H_idwt1d(Phi,z,wav,dwt_L,my,mx);
Ht=@(z) Ht_dwt1d(Phit,z,wav,dwt_L,my,mx);

W=@(z) midwt(z,wav,dwt_L);
Wt=@(z) mdwt(z,wav,dwt_L);

phantom_coeff=Wt(Img2D);
r_init=length(find(abs(phantom_coeff(:))>0));    %sparsity level for hard thresholding methods


% %solve by Back Projection
% s_BackProj=Ht(y);
% s_BackProj=reshape(s_BackProj,[my mx]);
% Img2D_BackProj=W(s_BackProj);
% Img2D_BackProj=Img2D_BackProj+Img2Dbar;
% PSNR_BackProj=psnr(Img2D,Img2D_BackProj,scale);
% 

%solve by DORE

z_init=zeros(m,1);
tic;
[s_DORE,A_index_DORE,Count_DORE]=DORE01(H,Ht,y,r_init,'Thresh',thresh,'visibility',0);
% 停止标准: 相对误差thresh, 最大迭代次数1000
t_DORE=toc
Count_DORE
s_DORE=reshape(s_DORE,[my mx]);
Img2D_DORE=W(s_DORE);
Img2D_DORE=Img2D_DORE+Img2Dbar;
PSNR_DORE=psnr(Img2D,Img2D_DORE,scale)

%solve by IDET
tic;
[s_IDET Out] = IDET(H, y, r_init, 'At', Ht, 'MaxIt',20,'Tolerance', thresh);
t_IDET = toc
Count_IDET = Out.iter
s_IDET=reshape(s_IDET,[my mx]);
Img2D_IDET=W(s_IDET);
Img2D_IDET=Img2D_IDET+Img2Dbar;
PSNR_IDET=psnr(Img2D,Img2D_IDET,scale)

%solve by IDETbeta
tic;
[s_IDETbeta Out] = IDETbeta(H, y, 0.8, 'At', Ht, 'MaxIt',20,'Tolerance', thresh);
% bata >0.8
t_IDETbeta = toc
Count_IDETbeta = Out.iter
s_IDETbeta=reshape(s_IDETbeta,[my mx]);
Img2D_IDETbeta=W(s_IDETbeta);
Img2D_IDETbeta=Img2D_IDETbeta+Img2Dbar;
PSNR_IDETbeta=psnr(Img2D,Img2D_IDETbeta,scale)

%solve by IDETgamma
tic;
[s_IDETgamma Out] = IDETgamma(H, y, 0.8, 'At', Ht, 'MaxIt',20,'Tolerance', thresh);
% gamma 0.5-0.8
t_IDETgamma = toc
Count_IDETgamma = Out.iter
s_IDETgamma=reshape(s_IDETgamma,[my mx]);
Img2D_IDETgamma=W(s_IDETgamma);
Img2D_IDETgamma=Img2D_IDETgamma+Img2Dbar;
PSNR_IDETgamma=psnr(Img2D,Img2D_IDETgamma,scale)

% %solve by SWGP
% tic;
% [s_SWGP Out] = SWGP(H, y, 0.8, 'At', Ht, 'MaxIt',20,'Tolerance', thresh);
% % gamma 0.5-0.8
% t_SWGP = toc
% Count_SWGP = Out.iter
% s_SWGP =reshape(s_SWGP,[my mx]);
% Img2D_SWGP =W(s_SWGP);
% Img2D_SWGP =Img2D_SWGP+Img2Dbar;
% PSNR_SWGP =psnr(Img2D,Img2D_SWGP,scale)


%Plotting
figure(1)
subplot(3,3,2)
imagesc(Img2Dcntr)
colormap(gray)
box off
axis off
axis equal
title('Original phantom image', 'fontsize',6);


subplot(3,3,3)
imagesc(smpl_pattern);
colormap(gray);
box off
axis off
axis equal
title('The 44 radial lines frequency sampling', 'fontsize',6);


% subplot(3,3,3)
% imagesc(Img2D_BackProj)
% colormap(gray)
% box off
% axis off
% title(['Back projectiom recovery (PSNR=',num2str(PSNR_BackProj),')'], 'fontsize',6);
% 
% 
% subplot(3,3,4)
% imagesc(Img2D_GPSR0)
% colormap(gray)
% box off
% axis off
% title(['(d) GPSR_0 recovery (PSNR=',num2str(PSNR_GPSR),')'], 'fontsize',6);
% 
% 
% subplot(3,3,5)
% imagesc(Img2D_GPSR)
% colormap(gray)
% box off
% axis off
% title(['(e) GPSR recovery (PSNR=',num2str(PSNR_GPSR),')'], 'fontsize',6);


% subplot(3,3,3)
% imagesc(Img2D_ADORE)
% colormap(gray)
% box off
% axis off
% title(['DORE recovery (PSNR=',num2str(PSNR_ADORE),')'], 'fontsize',6);


subplot(3,3,5)
imagesc(Img2D_DORE)
colormap(gray)
box off
axis off
axis equal
title(['DORE: PSNR=',num2str(PSNR_DORE,'%2.1f'), '; CPU time=',num2str(t_DORE,'%2.2f'), 's'], 'fontsize',6);

subplot(3,3,6)
imagesc(Img2D_IDET)
colormap(gray)
box off
axis off
axis equal
title(['IDET-k: PSNR=',num2str(PSNR_IDET,'%2.1f'), '; CPU time=',num2str(t_IDET,'%2.2f'), 's'], 'fontsize',6);

subplot(3,3,8)
imagesc(Img2D_IDETbeta)
colormap(gray)
box off
axis off
axis equal
title(['IDET-\beta=0.8: PSNR=',num2str(PSNR_IDETbeta,'%2.1f'), '; CPU time=',num2str(t_IDETbeta,'%2.2f'), 's'], 'fontsize',6);

subplot(3,3,9)
imagesc(Img2D_IDETgamma)
colormap(gray)
box off
axis off
axis equal
title(['IDET-\gamma=0.8: PSNR=',num2str(PSNR_IDETgamma,'%2.1f'), '; CPU time=',num2str(t_IDETgamma,'%2.2f'), 's'], 'fontsize',6);





if RadialNum == 44
    if imgsize == 128
        print -depsc2 Fig_phantom128
    else
        print -depsc2 Fig_phantom256
    end
end


PSNR  = [PSNR_DORE, PSNR_IDET, PSNR_IDETbeta, PSNR_IDETgamma];
Count = [Count_DORE, Count_IDET, Count_IDETbeta, Count_IDETgamma];
t     = [t_DORE, t_IDET, t_IDETbeta, t_IDETgamma];


save([img_name,num2str(my),'by',num2str(mx),'_RadialNum',num2str(RadialNum),'_thresh',num2str(thresh),'.mat']);

disp(sprintf('start at [%s]',mytime));
disp(sprintf('over at  [%s]',datestr(now)));

