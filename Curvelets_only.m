function [u_k,info] = Curvelets_only(image,nScales,nAngles,beta,lambda,maxIter,tol,sigma)
%TV-Curvelets Infimal Convolution.  
%INPUTS:
%image - '1' for phantom, '2' for medical image. 
%nScales - Number of scales.
%nAngles - Number of angles in the second coarse level.
%beta - The coefficient of the Curvelet term.
%lambda - The augmented Lagrangian term of the constraint v=Psi*u.
%maxIter - The maximum number of outer iterations.
%tol - for image '1' we use norm of relative true error, for image '2' we
%use norm of relative error from previous iteration.

clc, close all

%Info for output
info.normL2DataFit = [];
info.prim_res = [];
info.beta = beta;
info.lambda = lambda;
info.pcgIter = [];
info.rel_res_errors = [];
switch image
    case 1
        info.true_errors = [];
end
info.original_image = [];
info.noisy_image = [];
info.reconstruction = [];

%% Forming our image.

switch image
    case 1 %Phantom Image.
        Nx = 256;  % voxel number along x direction - outer square
        Ny = 256;  % voxel number along y direction - outer square
        nx = Nx/2; % voxel number along x direction - inner square
        ny = Ny/2; % voxel number along y direction - innerer square

        % specify a Curvelet in the inner square
        real = 1;      % real Curvelet transform
        finest = 1;    % all Curvelet

        % forward Curvelet transform, Psi: x -> y (vectorised)
        Psi = @(x) vectCurvelet(fdct_wrapping(x,real,finest,nScales,nAngles));            
        %C structure 
        C0 = fdct_wrapping(zeros(Nx,Ny),real,finest,nScales,nAngles); 
        C0vec = vectCurvelet(C0);
        Nc = length(C0vec);
        iPsi = @(y) ifdct_wrapping(unvectCurvelet(y,C0),real,Nx,Ny);
        Psi2 = @(x1,x2) [Psi(x1); Psi(x2)];
        iPsi2 = @(y1,y2) [iPsi(y1); iPsi(y2)];        
        plotCurvelet = @(c) fdct_wrapping_dispcoef(unvectCurvelet(c,C0));
        
        % specify the 24th Curvelet from top left in the 4th scale
        s = 2; %4 
        w = 24;
        C1 = C0;
        [A,B] = size(C1{s}{w});
        a = ceil((A+1)/2); % vertical location
        b = ceil((B+1)/2); % horizontal location
        C1{s}{w}(a,b) = sqrt(2)*10; % define the magnitude of the Curvelet in Fourier Domain
        C1vec = vectCurvelet(C1);
        disp(['||C1vec||_2 = ' num2str(norm(C1vec))])
        
        % backward Curvelet transform
        invC1 = iPsi(C1vec); %ifdct_wrapping(C1,real,Nx,Ny);
        disp(['||iPsi(C1vec)||_2 = ' num2str(norm(invC1(:)))])

        % construct phantom with single Curvelet in the inner square
        p0 = 1/2.*ones(Nx,Ny);
        p0(Nx/4+1:Nx/4+Nx/2,Ny/4+1:Ny/4+Ny/2) = zeros(nx,ny);
        p0 = p0 + invC1;

        % display the phantom
        figure(1); subplot(1,2,1);imagesc(p0); colorbar;title('Original Phantom')
        info.original_image = p0;
        
        % Construct noisy image.
        bottom = min(min(p0));
        top = max(max(p0)); 
        info.sigma = sigma;
        noise = sigma*randn(size(p0))*max(abs(p0(:)));
        p0_noise = p0 + noise;
        info.noisy_image = p0_noise;
        figure(1); subplot(1,2,2); imagesc(p0_noise); title("Noisy Image: \sigma = " + sigma); 
        caxis manual; caxis([bottom top]); colorbar;
        
    case 2 
        %Medical Image.
        p0 = im2double(imread('2DFilteredPalm.png'));
        p0 = imresize(p0,[256,256]);
        p0 = double(rgb2gray(imread('vesselTestImage.png')))/255;
        info.original_image = p0;
        % display the phantom
        figure(1); subplot(1,2,1);imagesc(p0); colorbar;title('Original Phantom')
        
        %Construct Noisy Image.
        bottom = min(min(p0));
        top = max(max(p0)); 
        info.sigma = sigma;
        noise = sigma*randn(size(p0))*max(abs(p0(:)));
        p0_noise = p0 + noise;
        info.noisy_image = p0_noise;
        figure(1); subplot(1,2,2); imagesc(p0_noise); title("Noisy Image: \sigma = " + sigma); 
        caxis manual; caxis([bottom top]); colorbar;
        size_p0 = size(p0);
        Nx = size_p0(1);
        Ny = size_p0(2);
        
        % specify a Curvelet in the inner square
        real = 1;      % real Curvelet transform
        finest = 1;    % all Curvelet

        % forward Curvelet transform, Psi: x -> y (vectorised)
        Psi = @(x) vectCurvelet(fdct_wrapping(x,real,finest,nScales,nAngles));            
        %C structure 
        C0 = fdct_wrapping(zeros(Nx,Ny),real,finest,nScales,nAngles); 
        C0vec = vectCurvelet(C0);
        Nc = length(C0vec);
        iPsi = @(y) ifdct_wrapping(unvectCurvelet(y,C0),real,Nx,Ny);
        Psi2 = @(x1,x2) [Psi(x1); Psi(x2)];
        iPsi2 = @(y1,y2) [iPsi(y1); iPsi(y2)];        
        plotCurvelet = @(c) fdct_wrapping_dispcoef(unvectCurvelet(c,C0));
end


%% Intial Vectors. 
f = p0_noise; %Matrix. Our starting Image.
u_k = f; %Matrix. Our solution Image.
v_length = length(Psi(u_k)); 
v = zeros(v_length,1); %Vector.
b = v; %Vector. 

%% Split Bregman Iteration.

bottom = min(min(p0));
top = max(max(p0));
figure(2); movegui(figure(2),'east');
subplot(2,3,1); imagesc(p0_noise); caxis manual; caxis([bottom top]); colorbar; title('Noisy Phantom');

n_iter = 0;

rel_res_error = 100;
while rel_res_error > tol && n_iter < maxIter
    
    n_iter = n_iter + 1;
    fprintf("OUTER Iteration = " + n_iter + "\n");
    
    %u-update, using differentiation.
    u_k_1 = (f + lambda*iPsi(b + v))/(1+lambda);
    info.normL2DataFit(n_iter) = norm(f(:) - u_k_1(:),2);
    
    rel_res_error = norm(u_k_1 - u_k)/norm(u_k_1)
    info.rel_res_errors(n_iter) = rel_res_error;
    true_error = norm(u_k_1 - p0);
    info.true_errors(n_iter) = true_error;


    figure(2); subplot(2,3,2); imagesc(u_k_1); caxis manual; caxis([bottom top]); colorbar; 
    title("Reconstructed Image,Iteration=" + n_iter + "\newline \lambda=" ...
        + lambda + ",\beta=" + beta); 
    
    %v-update.
    d = b - Psi(u_k_1);
    v = softThresh(-d,beta/lambda);
    iCv = iPsi(v); 
    figure(2); 
    subplot(2,3,5); imagesc(iCv); colorbar; title("iCv, Iteration = " + n_iter);drawnow;
    subplot(2,3,6); imagesc(abs(fftshift(fft2(ifftshift(iCv))))); colorbar; title("Cv, Iteration = " + n_iter);drawnow;
    info.normL1v(n_iter) = norm(v(:),1);
    
    %Primal-Residual Update;
    info.prim_res(n_iter) = norm(v - Psi(u_k_1));
    
    %b-update.
    b = b + v - Psi(u_k_1);
    figure(3); movegui(figure(3),'east');
    subplot(1,3,3); imagesc(plotCurvelet(b)); colorbar; title("b, Iteration = " + n_iter); drawnow;


    %Update iterations.
    u_k = u_k_1;
    
    %Collect info about iteration
    info.us(:,n_iter) = u_k(:);
    info.vs(:,n_iter) = v(:);
    info.objective(n_iter) = 0.5*info.normL2DataFit(n_iter) + beta*info.normL1v(n_iter);

    
end

%% Final Results.

info.reconstruction = u_k_1;

switch image
    case 1
        figure(4); 
        subplot(2,2,1); imagesc(p0); axis image; caxis manual; caxis([bottom top]); colorbar; 
        title('Original Phantom')
        subplot(2,2,2); imagesc(p0_noise); axis image; caxis manual; caxis([bottom top]); colorbar; 
        title('Noisy Phantom')
        subplot(2,2,3); imagesc(u_k_1);  axis image; 
        title("Reconstructed Image,Iteration=" + n_iter + "\newline \lambda=" ...
        + lambda + ",\beta=" + beta); 
        caxis manual; caxis([bottom top]); colorbar;
        subplot(2,2,4); imagesc(p0 - u_k_1);  axis image; 
        title("Difference between \newline reconstruction and original image"); colorbar
        figure; semilogy(info.true_errors); title('True Error');
        xlabel('Iteration');
        figure; semilogy(info.prim_res); title('Primal Residual');
        xlabel('Iteration');
    case 2
        figure; 
        imagesc(p0_noise); axis image; caxis manual; caxis([bottom top]);
        title("Noisy Image, \sigma=" + sigma); colorbar;
        
        figure; 
        imagesc(p0); axis image; caxis manual; caxis([bottom top]);
        title("Original Image"); colorbar;
        figure;
        imagesc(u_k_1); axis image; caxis manual; caxis([bottom top]);
        title("Final \lambda=" ...
        + lambda + ",\beta=" + beta + ...
        ",noise=" + sigma); 
        colorbar;
        figure; semilogy(info.true_errors); title('True Error');
        xlabel('Iteration');
        figure; semilogy(info.prim_res); title('Primal Residual');
        xlabel('Iteration');
end

figure;
semilogy(info.objective, '-b'); hold on; 
semilogy(0.5*info.normL2DataFit, '.-k'); 
semilogy(info.beta*info.normL1v, '--g');
title("Objective function"); legend('Objective(k)', '$\|f-u_k\|_2$', '$\beta \|v_k\|_1$', 'Interpreter','Latex');
hold off



%% function handles for vectorizing/unvectorizing the Curvelet coefficients
function x = vectCurvelet(C)
% VECTCURVELET Vectorizes Curvelet coefficients, ordering: scale, angle, (:)
%  x = vectCurvelet(C)

i = 0;
for s = 1:length(C) %loop through scales
  for w = 1:length(C{s}) %loop through angles
    lCsw = length(C{s}{w}(:));
    x(i + (1:lCsw)) = C{s}{w}(:);
    i = i + lCsw;
  end
end
x = x'; %return column vector
end


function C = unvectCurvelet(x, S)
% UNVECTCURVELET Reshapes vector of Curvelet coefficients x to Curvelet with structure S
%  C = unvectCurvelet(x, S)

C = S; %initialize with S for speed

i = 0;
for s = 1:length(S) %loop through scales
  for w = 1:length(S{s}) %loop through angles
    sizeCsw = size(S{s}{w});
    C{s}{w} = reshape(x(i + (1:prod(sizeCsw))), sizeCsw(1), sizeCsw(2));
    i = i + prod(sizeCsw);
  end
end
end


end



