function [u_k,info] = TVCurvelets_21_8(image,u_solver,nScales,nAngles,alpha,beta,lambda,mu,maxIter,non_negativity,tol,sigma)
%TV-Curvelets Infimal Convolution.  
%INPUTS:
%image - '1' for phantom, '2' for medical image. 
%u_solver - '1' for pcg solver. '2' for discrete cosine matrix solver.
%nScales - Number of scales.
%nAngles - Number of angles in the second coarse level.
%alpha - The coefficient of the TV term.
%beta - The coefficient of the Curvelet term.
%lambda - The augmented Lagrangian term of the first constraint: A1u + A2z
%+ A3w + a4v = 0.
%mu - The augmented Lagrangian term of the second constraint: w=s.
%maxIter - The maximum number of outer iterations.
%non_negativity - '1' for no prior, '2' for non-negativity prior.
%tol - for image '1' we use norm of relative true error, for image '2' we
%use norm of relative error from previous iteration.
%sigma - percentage of noise of maximum absolute value.

clc, close all

%Info for output
info.normL2DataFit = [];
info.normL1z = [];
info.normL1w = [];
info.prim_res = [];
info.alpha = alpha;
info.beta = beta;
info.lambda = lambda;
info.mu = mu;
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
        p0_nc = p0;
        p0 = p0 + invC1;

        % display the phantom
        figure(1); subplot(1,2,1);imagesc(p0); colorbar;title('Original Phantom')
        info.original_image = p0;
        
        % forward Curvelet transform
        Cp0 = fdct_wrapping(p0,real,finest,nScales,nAngles);
        Cp0_nc = fdct_wrapping(p0_nc,real,finest,nScales,nAngles);

        % display the Curvelet coefficients
        img = fdct_wrapping_dispcoef(Cp0);
        figure;imagesc(abs(img));axis image;title('Curvelet Coeffs of p0 with curvelet')
        img = fdct_wrapping_dispcoef(Cp0_nc);
        figure;imagesc(abs(img));axis image;title('Curvelet Coeffs of p0 without curvelet')

        % Construct noisy image.
        bottom = min(min(p0));
        top = max(max(p0));
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
z1 = u_k; %Matrix.
z2 = u_k; %Matrix.
w1 = u_k; %Matrix.
w2 = u_k; %Matrix.
b1_1 = u_k; %Matrix.
b1_2 = u_k; %Matrix.
s = [u_k;u_k]; %Matrix.
a = s; %Matrix. 
v_length = length(Psi2(w1,w2)); 
v = zeros(v_length,1); %Vector.
b2 = v; %Vector. 

%% Split Bregman Iteration.

bottom = min(min(p0));
top = max(max(p0));
figure(2); movegui(figure(2),'east');
subplot(2,3,1); imagesc(p0_noise); caxis manual; caxis([bottom top]); colorbar; title('Noisy Phantom');

n_iter = 0;
[Dx,Dy] = LinOpTV(u_k(:));
D2 = Dx'*Dx + Dy'*Dy;
A1 = speye(size(D2)) + lambda*(D2);

switch u_solver
    case 1 
    case 2
        Ae1 = A1(:,1);
        e1 = zeros(length(D2),1);
        e1(1) = 1;
        D_top = dct2(full(Ae1));
        D_bottom = dct2(e1);
        A1_diag_vec = D_top./D_bottom;
        A1_diag_inv_vec = 1./A1_diag_vec;
        A1_diag_inv_mat = spdiags(A1_diag_inv_vec,0,Nx^2,Nx^2);
end

rel_res_error = 100;
while rel_res_error > tol && n_iter < maxIter
    
    n_iter = n_iter + 1;
    fprintf("OUTER Iteration = " + n_iter + "\n");
    
    %u-update.
    DxT_z_w_b = Dx'*(z1(:)+w1(:)+b1_1(:));
    DyT_z_w_b = Dy'*(z2(:)+w2(:)+b1_2(:));
    c1 = f(:) + lambda*(DxT_z_w_b + DyT_z_w_b);
    
    switch u_solver
        case 1
            fprintf('u PCG \n');
            [u_k_1,flag,relres,iter,resvec] = pcg(A1,c1,1e-6,100);
            info.pcgIter(n_iter) = iter;
            u_k_1 =  reshape(u_k_1,[Nx,Ny]);
        case 2
            term1 = A1_diag_inv_mat*dct2(c1);
            u_k_1 = idct2(term1);
            u_k_1 =  reshape(u_k_1,[Nx,Ny]);
    end
    
    rel_res_error = norm(u_k_1 - u_k)/norm(u_k_1)
    info.rel_res_errors(n_iter) = rel_res_error; 
    true_error = norm(u_k_1 - p0);
    info.true_errors(n_iter) = true_error; 

    figure(2); subplot(2,3,2); imagesc(u_k_1); caxis manual; caxis([bottom top]); colorbar; 
    title("Reconstructed Image,Iteration=" + n_iter + "\newline \lambda=" ...
        + lambda + ",\mu=" + mu + ",\alpha=" + alpha + ",\beta=" + beta); 
    
    %New vectors with new 'u'.
    Dxu_k_1_vec = Dx*(u_k_1(:)); 
    Dxu_k_1 = reshape(Dxu_k_1_vec,[Nx,Ny]); %Matrix
    Dyu_k_1_vec = Dy*(u_k_1(:)); 
    Dyu_k_1 = reshape(Dyu_k_1_vec,[Nx,Ny]); %Matrix
    Du_k_1_vec = [Dxu_k_1_vec;Dyu_k_1_vec];
    Du_k_1 = [Dxu_k_1;Dyu_k_1]; %Matrix
    info.normL2DataFit(n_iter) = norm(f(:) - u_k_1(:),2);
    
    %z-update.
    d1_z = -Du_k_1_vec + [w1(:);w2(:)] + [b1_1(:);b1_2(:)];
    z_vec = softThresh(-d1_z,alpha/lambda);
    z1 = reshape(z_vec(1:Nx*Ny),[Nx,Ny]);
    z2 = reshape(z_vec((Nx*Ny)+1:end),[Nx,Ny]);
    figure(2); subplot(2,3,3); imagesc([z1;z2]); colorbar; title("[z1;z2], Iteration = " + n_iter);drawnow;
    info.normL1z(n_iter) = norm(z_vec(:),1);
    
    %w-update, using differentiation.
    d1_w = -Du_k_1 + [z1;z2] + [b1_1;b1_2];
    d2_w = v + b2;
    switch non_negativity
        case 1
            w = 0.5*(iPsi2(d2_w(1:Nc),d2_w(Nc+(1:Nc))) - d1_w);
        case 2
            w = (mu*(s - a) + lambda*(iPsi2(d2_w(1:Nc),d2_w(Nc+(1:Nc))) - d1_w))/(mu + 2*lambda);

    end
    w1 = w(1:Nx,1:end);
    w2 = w(Nx+1:end,1:end);
    figure(2);
    subplot(2,3,4); imagesc([w1;w2]); colorbar; title("[w1;w2], Iteration = " + n_iter);drawnow;
    
    %s-update.
    s = max(0,[w1;w2] + a);
    
    %a-update.
    a = a + [w1;w2] - s;
    
    %v-update.
    d2_v = -Psi2(w1, w2) + b2;
    v = softThresh(-d2_v,beta/lambda);
    iCv = iPsi2(v(1:Nc), v(Nc+(1:Nc))); 
    iCv1 = iCv(1:Nx,1:end);
    iCv2 = iCv(Nx+1:end,1:end);
    figure(2); movegui(figure(2),'east');
    subplot(2,3,5); imagesc([iCv1;iCv2]); colorbar; title("[iCv1;iCv2], Iteration = " + n_iter);drawnow;
    subplot(2,3,6); imagesc([abs(fftshift(fft2(ifftshift(iCv1)))); abs(fftshift(fft2(ifftshift(iCv2))))]); colorbar; title("[Cv1;Cv2], Iteration = " + n_iter);drawnow;
    info.normL1v(n_iter) = norm(v(:),1);
    
    %Primal-Residual update.
    info.prim_res(n_iter) = norm([- Dxu_k_1 + z1 + w1;...
                            - Dyu_k_1 + z2 + w2])...
                            + norm(- Psi2(w1,w2) + v);
    
    %b-update.
    b1_1 = b1_1 - Dxu_k_1 + z1 + w1;
    b1_2 = b1_2 - Dyu_k_1 + z2 + w2;
    b2 = b2 - Psi2(w1,w2) + v;
    figure(3); 
    subplot(1,3,1); imagesc(b1_1); colorbar; title("[b1_1], Iteration = " + n_iter); drawnow;
    subplot(1,3,2); imagesc(b1_2); colorbar; title("[b1_2], Iteration = " + n_iter); drawnow;
    subplot(1,3,3); imagesc([plotCurvelet(b2(1:Nc)); plotCurvelet(b2(Nc+(1:Nc)))]); colorbar; title("[b2], Iteration = " + n_iter); drawnow;


    %Update iterations.
    u_k = u_k_1;
    
    %Collect info about iteration
    info.us(:,n_iter) = u_k(:);
    info.vs(:,n_iter) = v(:);
    info.ws(:,n_iter) = [w1(:);w2(:)];
    info.objective(n_iter) = 0.5*info.normL2DataFit(n_iter) + alpha*info.normL1z(n_iter) + beta*info.normL1v(n_iter);

    
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
        + lambda + ",\mu=" + mu + ",\alpha=" + alpha + ",\beta=" + beta); 
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
        + lambda + ",\mu=" + mu + ",\alpha=" + alpha + ",\beta=" + beta + ...
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
semilogy(info.alpha*info.normL1z, '-.r'); 
semilogy(info.beta*info.normL1v, '--g');
title("Objective function"); legend('Objective(k)', '$\|f-u_k\|_2$', '$\alpha \|z_k\|_1$', '$\beta \|v_k\|_1$', 'Interpreter','Latex');
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



