%% Inverse_problem

% getting x from Ay

%% Diffuse optics Green's function for an infinite medium, computed in a rectangular slab
load('NeuroDOT_Data_Sample_CCW1.mat'); % data, info, flags


mua = .0192; % flags.op.mua_gray=[0.0180,0.0192];
musp = 0.6726; % flags.op.musp_gray=[0.8359,0.6726];
nu = 1.4; 
D = 1/(3*(mua+musp)); % diffusion coefficient
mu_eff = sqrt(mua/D); % effective attenuation coefficien

% Define bounds on medium
xBnds = [-20 20]; yBnds = [-30 30]; zBnds = [1 30];  
mmX = 2; mmY = 2; mmZ = 2; 

[Y X Z] = meshgrid(yBnds(1):mmY:yBnds(2), xBnds(1):mmX:xBnds(2), zBnds(1):mmZ:zBnds(2)); % generate coordinates for slab 

voxCrd = double([X(:) Y(:) Z(:)]); % coordinates (mm) for each voxel, reshaped as 1D vector

nX = size(X,1); nY = size(X,2); nZ = size(X,3); % volume is nX x nY x nZ voxels
srcPos =info.optodes.spos3;
detPos =info.optodes.dpos3;


%% A Matrix calculation
r = pdist2(srcPos, voxCrd); % distance from source to each voxel
GsAnalytic = 1./(4*pi*D*r).*exp(-mu_eff*r); % Green's function, source to voxel

rdet = pdist2(voxCrd, detPos);
GsDet = 1./(4*pi*D*rdet).*exp(-mu_eff*rdet); % Green's function, voxel to detector

r2 = pdist2(detPos,srcPos );
Gs = 1./(4*pi*D*r2).*exp(-mu_eff*r2);

nsrc =length(srcPos);
ndet =length(detPos);
A = [];
A_tmp =[];
for i = 1:nsrc
    for j = 1: ndet

        A_tmp = GsAnalytic(i,:).* GsDet(:,j)';
        A = [A; A_tmp];
    end
end

Gs = reshape(Gs, size(Gs, 1) * size(Gs, 2), []);
A = A ./ Gs;




%% temporal evolution though time
freq = 1; %hz
x = zeros(size(A,2),1);
x(1) = 1;
y = A*x;
lambda = 0.0001; %max(A,[],"all")/400;

A_reg = A'*(A*A'+lambda*eye(size(A,1)))^-1;
% A_reg = (A'*A+lambda*eye(size(A,1)))\A';
x_est = A_reg *y;

figure;
imagesc(A_reg*A);
%%
reconstrution = reshape(x_est,nX,nY,nZ);

figure;
plot(x');
hold on
plot(x_est');
ylabel("intentsity changes");
xlabel("tissue structures");
legend("x","x_estimate");

% figure;
% imagesc(squeeze(reconstrution(:,:,15)));
% colorbar;


reconstrutionOriginal = reshape(x,nX,nY,nZ);

figure;
sliceViewer(reconstrutionOriginal);

figure;
sliceViewer(reconstrution);
colorbar;

figure;
% plot(1:length(y), y);
imagesc(reconstrutionOriginal(:, :, 1));