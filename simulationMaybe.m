%% Diffuse optics Green's function for an infinite medium, computed in a rectangular slab
load('NeuroDOT_Data_Sample_CCW1.mat'); % data, info, flags

mua = .0192; % flags.op.mua_gray=[0.0180,0.0192];
musp = 0.6726; % flags.op.musp_gray=[0.8359,0.6726];
nu = 1.4; 
D = 1/(3*(mua+musp));
mu_eff = sqrt(mua/D); 

% Define bounds on medium
% xBnds = [-30 30]; yBnds = [-45 45]; zBnds = [1 30];  
xBnds = [-75 75]; yBnds = [-45 45]; zBnds = [1 80];  
mmX = 2; mmY = 2; mmZ = 2; 

[Y X Z] = meshgrid(yBnds(1):mmY:yBnds(2), xBnds(1):mmX:xBnds(2), zBnds(1):mmZ:zBnds(2)); % generate coordinates for slab 

voxCrd = double([X(:) Y(:) Z(:)]); % coordinates (mm) for each voxel, reshaped as 1D vector

nX = size(X,1); nY = size(X,2); nZ = size(X,3); % volume is nX x nY x nZ voxels

srcPos = [0 0 0; 0 0 1]; % source position, in 3D coordinates 
detPos = [2 2 2; 3 3 3; 4 4 4];
detPos = [info.optodes.dpos3(:,1),info.optodes.dpos3(:,2),info.optodes.dpos3(:,3)];
srcPos = [info.optodes.spos3(:,1),info.optodes.spos3(:,2),info.optodes.spos3(:,3)];


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

regions = 10;
s = floor(length(perturbations)/regions);

perturbations = zeros(size(voxCrd, 1), 1);
measurements = zeros(size(A, 1), regions);

for i = 1:regions
    perturbations = zeros(size(voxCrd, 1), 1);
    perturbations(s * i : min(s * i + s - 1, length(perturbations))) = 1;
    measurements(:, i) = A * perturbations;
end

figure, imagesc(log(measurements));
% A = GsAnalytic * GsDet;

% # measurements x number of voxels, ie 2 detector, 3 detectors, 6
% measurements, ie 6 x 8Meg

% tmp = reshape(GsAnalytic,nX,nY,nZ); % reshape as 3D volume
% figure, sliceViewer(tmp,'Colormap',hot(256)); % simple viewer
% 
% figure, sliceViewer(log10(tmp),'Colormap',hot(256)); % simple viewer, log compressed
% 
