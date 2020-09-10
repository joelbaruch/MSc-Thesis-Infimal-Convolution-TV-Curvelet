%%
clc, close all, clear all
tic
%TVCurvelets(image,u_solver,nScales,nAngles,alpha,beta,lambda,mu,maxIter,non_negativity,tol,sigma)
[u_k, info] = TVCurvelets_21_8(1,1,4,128,0.1,0.1,10,20,100,1,1e-5,0.1);
toc
alpha = info.alpha;
beta = info.beta;
psnr_recon = psnr(info.reconstruction,info.original_image);
psnr_noisy = psnr(info.noisy_image,info.original_image);
ssim_recon = ssim(info.reconstruction,info.original_image);
ssim_noisy = ssim(info.noisy_image,info.original_image); 
%%
clc, close all, clear all
tic
%TV_only(image,u_solver,nScales,nAngles,alpha,lambda,maxIter,tol,sigma)
[u_k,info] = TV_only(2,1,4,128,0.02,10,100,1e-5,0.03);
toc
alpha = info.alpha;
psnr_recon = psnr(info.reconstruction,info.original_image);
psnr_noisy = psnr(info.noisy_image,info.original_image);
ssim_recon = ssim(info.reconstruction,info.original_image);
ssim_noisy = ssim(info.noisy_image,info.original_image); 
%%
clc, close all, clear all
tic
%Curvelets_only(image,nScales,nAngles,beta,lambda,maxIter,tol,sigma)
[u_k,info] = Curvelets_only(2,4,128,0.01,10,100,1e-5,0.03);
toc
beta = info.beta;
psnr_recon = psnr(info.reconstruction,info.original_image);
psnr_noisy = psnr(info.noisy_image,info.original_image);
ssim_recon = ssim(info.reconstruction,info.original_image);
ssim_noisy = ssim(info.noisy_image,info.original_image); 