%% Diffuse optics Green's function for an infinite medium, computed in a rectangular slab

mua = .0192; % flags.op.mua_gray=[0.0180,0.0192];
musp = 0.6726; % flags.op.musp_gray=[0.8359,0.6726];
nu = 1.4; 
D = 1/(3*(mua+musp)); % diffusion coefficient
mu_eff = sqrt(mua/D); % effective attenuation coefficien

% Define bounds on medium
xBnds = [-30 30]; yBnds = [-45 45]; zBnds = [1 30];  
mmX = 2; mmY = 2; mmZ = 2; 

[Y X Z] = meshgrid(yBnds(1):mmY:yBnds(2), xBnds(1):mmX:xBnds(2), zBnds(1):mmZ:zBnds(2)); % generate coordinates for slab 

voxCrd = double([X(:) Y(:) Z(:)]); % coordinates (mm) for each voxel, reshaped as 1D vector

nX = size(X,1); nY = size(X,2); nZ = size(X,3); % volume is nX x nY x nZ voxels

srcPos = [0 0 0]; % source position, in 3D coordinates 

r = pdist2(srcPos,voxCrd); % distance from source to each voxel

GsAnalytic = 1./(4*pi*D*r).*exp(-mu_eff*r); % Green's function, each voxel

tmp = reshape(GsAnalytic,nX,nY,nZ); % reshape as 3D volume

% figure, sliceViewer(tmp,'Colormap',hot(256)); % simple viewer

figure, sliceViewer(log10(tmp),'Colormap',hot(256)); % simple viewer, log compressed
%% Detection source function

% GsAnalytic = green_function(voxCrd,srcPos); %detext corrd
srcPos= [-20 -30 10; 20 30 10];
detPos= [0,0,0];

dis_btw_src = pdist2(srcPos(1,:),srcPos(2,:));
disp(dis_btw_src);

r = pdist2(srcPos, voxCrd); % distance from source to each voxel
GsAnalytic = 1./(4*pi*D*r).*exp(-mu_eff*r); % Green's function, source to voxel

rdet = pdist2(voxCrd, detPos);
GsDet = 1./(4*pi*D*rdet).*exp(-mu_eff*rdet); % Green's function, voxel to detector

r2 = pdist2(srcPos, detPos);
Gs = 1./(4*pi*D*r2).*exp(-mu_eff*r2);

GsAnalytic = reshape(GsAnalytic, size(GsAnalytic, 1), 1, []);
GsDet = reshape(GsDet, 1, size(GsDet, 2), []);

A = GsAnalytic .* GsDet;
A = reshape(A, size(A, 1) * size(A, 2), []);
Gs = reshape(Gs, size(Gs, 1) * size(Gs, 2), []);
A = A ./ Gs;
A = A(1,:)+A(2,:);

% show A matrix in one voxel
% change the source detector pairs

%  brain activities: change in orgniization of area: changes signales
%  green functions indicate the change in absorbtion relative to background
%  Coupleing coeffeicent: vairance of the sources 

tmp = reshape(A,nX,nY,nZ); % reshape as 3D volume

figure, sliceViewer(log10(tmp),'Colormap',hot(256)); % simple viewer, lo
title(sprintf('Source Positions: [%d, %d, %d; %d, %d, %d]', ...
    srcPos(1,1), srcPos(1,2), srcPos(1,3), srcPos(2,1), srcPos(2,2), srcPos(2,3)));
colorbar;
%% functions

function GsAnalytic = green_function(pos_1, pos_2)
    mua = .0192; % flags.op.mua_gray=[0.0180,0.0192];
    musp = 0.6726; % flags.op.musp_gray=[0.8359,0.6726];
    nu = 1.4; 
    D = 1/(3*(mua+musp)); % diffusion coefficient
    mu_eff = sqrt(mua/D); % effective attenuation coefficien

    r = pdist2(pos_1, pos_2); % distance from source to each voxel
    GsAnalytic = 1./(4*pi*D*r).*exp(-mu_eff*r); % Green's function, each voxel

    % nX = dim(1); nY = dim(2); nZ = dim(3); % volume is nX x nY x nZ voxels
    % tmp = reshape(GsAnalytic,nX,nY,nZ); % reshape as 3D volume

end


