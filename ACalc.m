%% Diffuse optics Green's function for an infinite medium, computed in a rectangular slab

mua = .0192; % flags.op.mua_gray=[0.0180,0.0192];
musp = 0.6726; % flags.op.musp_gray=[0.8359,0.6726];
nu = 1.4; 
D = 1/(3*(mua+musp));
mu_eff = sqrt(mua/D); 

% Define bounds on medium
xBnds = [-30 30]; yBnds = [-45 45]; zBnds = [1 30];  
mmX = 2; mmY = 2; mmZ = 2; 

[Y X Z] = meshgrid(yBnds(1):mmY:yBnds(2), xBnds(1):mmX:xBnds(2), zBnds(1):mmZ:zBnds(2)); % generate coordinates for slab 

voxCrd = double([X(:) Y(:) Z(:)]); % coordinates (mm) for each voxel, reshaped as 1D vector

nX = size(X,1); nY = size(X,2); nZ = size(X,3); % volume is nX x nY x nZ voxels

srcPos = [0 0 0; 0 0 1]; % source position, in 3D coordinates 
detPos = [2 2 2; 3 3 3; 4 4 4];


r = pdist2(srcPos,voxCrd); % distance from source to each voxel

GsAnalytic = 1./(4*pi*D*r).*exp(-mu_eff*r); % Green's function, each voxel

rdet = pdist2(voxCrd, detPos);
GsDet = 1./(4*pi*D*rdet).*exp(-mu_eff*rdet); % Green's function, each voxel

r2 = pdist2(srcPos, detPos);
Gs = 1./(4*pi*D*r2).*exp(-mu_eff*r2);

GsAnalytic = reshape(GsAnalytic, size(GsAnalytic, 1), 1, []);
GsDet = reshape(GsDet, 1, size(GsDet, 2), []);

A = GsAnalytic .* GsDet;
A = reshape(A, size(A, 1) * size(A, 2), []);
Gs = reshape(Gs, size(Gs, 1) * size(Gs, 2), []);
A = A ./ Gs;

% A = GsAnalytic * GsDet;

% # measurements x number of voxels, ie 2 detector, 3 detectors, 6
% measurements, ie 6 x 8Meg

% tmp = reshape(GsAnalytic,nX,nY,nZ); % reshape as 3D volume
% figure, sliceViewer(tmp,'Colormap',hot(256)); % simple viewer
% 
% figure, sliceViewer(log10(tmp),'Colormap',hot(256)); % simple viewer, log compressed
% 
