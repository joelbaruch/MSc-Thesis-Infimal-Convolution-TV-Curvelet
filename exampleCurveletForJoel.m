% THIS FILE IS AN EXAMPLE ONLY. DO NOT MODIFY THIS FILE TO RUN EXPERIMENTS! 
% MAKE YOUR OWN COPY IN WHICH YOU ADJUST PATHS, OPTIONS ETC.
%
% Copyright (C) 2020 Bolin Pan

clc; close all; clear all



%% initialization
Nx = 256;  % voxel number along x direction - outer square
Ny = 256;  % voxel number along y direction - outer square
nx = Nx/2; % voxel number along x direction - inner square
ny = Ny/2; % voxel number along y direction - innerer square

% p0Mask = ones(Nx,Ny); % outer square
% p0Mask(Nx/4+1:Nx/4+Nx/2,Ny/4+1:Ny/4+Ny/2) = 0; % inner square
% p0Mask = 1/2.*p0Mask; % defined the contrast


%% specify a Curvelet in the inner square
nScales = 4;   % number of scale
nAngles = 128; % number of angles in the 2nd coarse level
real = 1;      % real Curvelet transform
finest = 1;    % all Curvelet

% forward Curvelet transform
innerSquare = zeros(nx,ny);
C = fdct_wrapping(innerSquare,real,finest,nScales,nAngles);

% specify the 24th Curvelet from top left in the 4th scale
s = 4;
w = 24;
[A,B] = size(C{s}{w});
a = ceil((A+1)/2); % vertical location
b = ceil((B+1)/2); % horizontal location
C{s}{w}(a,b) = 10; % define the magnitude of the Curvelet in Fourier Domain

% backward Curvelet transform
innerSquareC = ifdct_wrapping(C,real,nx,ny);


%% construct phantom with single Curvelet in the inner square
p0 = 1/2.*ones(Nx,Ny);
p0(Nx/4+1:Nx/4+Nx/2,Ny/4+1:Ny/4+Ny/2) = innerSquareC;

% display the phantom
figure;imagesc(p0);axis image;colorbar;title('phantom')


%% forward backward transform on the phantom
nScalesp0 = 4;   % default would be ceil(log2(min(Nx,Ny)) - 3)]
nAnglesp0 = 128; % default minimum 8
real = 1;      % real Curvelet transform
finest = 1;    % all Curvelet

% forward Curvelet transform
Cp0 = fdct_wrapping(p0,real,finest,nScalesp0,nAnglesp0);
vCp0 = vectCurvelet(Cp0); % vectorize the coefficients

% display the Curvelet coefficients
img = fdct_wrapping_dispcoef(Cp0);
figure;imagesc(abs(img));axis image;title('Curvelet Coeffs of p0')

% backward Curvelet transform
mask = fdct_wrapping(ones(Nx,Ny),real,finest,nScalesp0,nAnglesp0); % construct mask
uvCp0 = unvectCurvelet(vCp0, mask); % unvectorized the coefficients
p0Recon = ifdct_wrapping(uvCp0,real,Nx,Ny);

% display the reconstruction
figure;imagesc(p0Recon);axis image;colorbar;title('recovered phantom')
figure;imagesc(p0Recon-p0);axis image;colorbar;title('recovered error')



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




