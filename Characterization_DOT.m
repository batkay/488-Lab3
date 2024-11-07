close all
%% light fall-offs
subplot(1,2,1)
semilogy(data(info.pairs.r3d < 30,:)')
title("Data Traces, source-detector distance < 30 mm ")
ylabel('Intensity, \muW'), xlabel('t, samples')

subplot(1,2,2)
semilogy(info.pairs.r3d,mean(data,2),'.')
title("Light Falloff")
xlabel("source-detector distance, mm")
ylabel("Mean Signal Intensity")
%% Diffuse optics Green's function for an infinite medium, computed in a rectangular slab

mua = .0192; % flags.op.mua_gray=[0.0180,0.0192];
musp = 0.6726; % flags.op.musp_gray=[0.8359,0.6726];
nu = 1.4; 
D = 1/(3*(mua+musp)); % diffusion coefficient
mu_eff = sqrt(mua/D); % effective attenuation coefficien

% Define bounds on medium
xBnds = [-60 60]; yBnds = [-80 80]; zBnds = [1 30];  
mmX = 2; mmY = 2; mmZ = 2; 

[Y X Z] = meshgrid(yBnds(1):mmY:yBnds(2), xBnds(1):mmX:xBnds(2), zBnds(1):mmZ:zBnds(2)); % generate coordinates for slab 

voxCrd = double([X(:) Y(:) Z(:)]); % coordinates (mm) for each voxel, reshaped as 1D vector

nX = size(X,1); nY = size(X,2); nZ = size(X,3); % volume is nX x nY x nZ voxels

%% Parameters for source and detector layout
numSources = 10;
numDetectors = 20;
sourceSpacing = 2; % distance between sources in mm
detectorSpacing = 5; % distance between detectors in mm


% placing sources along the x-axis
srcPos = [ (1:numSources)' * sourceSpacing, zeros(numSources, 1), zeros(numSources, 1) ];

% placing detectors along the y-axis at a set distance from the x-axis
detPos = [ zeros(numDetectors, 1), (1:numDetectors)' * detectorSpacing, zeros(numDetectors, 1) ];

% % Define the bounds for source and detector positions
% sourceBounds = [-30, 30]; % Define the range for sources in each dimension (mm)
% detectorBounds = [-50, 50]; % Define the range for detectors in each dimension (mm)
% 
% % Generate 10 random source positions within the source bounds
% numSources = 10;
% srcPos = sourceBounds(1) + (sourceBounds(2) - sourceBounds(1)) * rand(numSources, 3);
% 
% % Generate 20 random detector positions within the detector bounds
% numDetectors = 20;
% detPos = detectorBounds(1) + (detectorBounds(2) - detectorBounds(1)) * rand(numDetectors, 3);
% 

%% Parameters for source and detector layout on a curved surface
% numSources = 24;
% numDetectors = 28;
% radius = 100; % Radius of curvature in mm (approximate head radius)
% 
% % Define angular coverage
% angularCoverage = pi/2; % 90 degrees coverage
% 
% % Calculate angular positions for sources and detectors
% thetaSources = linspace(-angularCoverage/2, angularCoverage/2, numSources);
% thetaDetectors = linspace(-angularCoverage/2, angularCoverage/2, numDetectors);
% 
% % Calculate positions of sources and detectors on the arc
% srcPos = [radius * cos(thetaSources)', radius * sin(thetaSources)', zeros(numSources, 1)];
% detPos = [radius * cos(thetaDetectors)', -radius * sin(thetaDetectors)', zeros(numDetectors, 1)];
% 
% % Plot to visualize the layout
% figure;
% plot3(srcPos(:,1), srcPos(:,2),srcPos(:,3), 'ro', 'DisplayName', 'Sources'); hold on;
% plot3(detPos(:,1), detPos(:,2), detPos(:,3), 'bo', 'DisplayName', 'Detectors');
% legend show;
% axis equal;
% title('24x28 Source and Detector Positions on a Curved Surface');
% xlabel('X Position (mm)');
% ylabel('Y Position (mm)');
srcPos =info.optodes.spos3;
detPos =info.optodes.dpos3;


%% 

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
%% 

A_sum = sum(A,1);

tmp = reshape(A_sum,nX,nY,nZ); % reshape as 3D volume

figure, sliceViewer(log10(tmp),'Colormap',hot(256)); % simple viewer, lo
title(sprintf('Source Positions: [%d, %d, %d; %d, %d, %d]', ...
    srcPos(1,1), srcPos(1,2), srcPos(1,3), srcPos(2,1), srcPos(2,2), srcPos(2,3)));
colorbar;

%% 
perturbations = zeros(size(voxCrd, 1), 1);
perturbations(1) = 1;

measurements = A * perturbations;

% Compute the mean intensity for each source-detector pair in A
meanIntensity = mean(A, 2); % mean across all voxels for each source-detector pair

% Extract source-detector distances from r2
sourceDetectorDistances = r2(:); % flatten the matrix to a vector

% % Plot the light fall-off
% figure;
% semilogy(sourceDetectorDistances, meanIntensity, 'o-'); % use 'o-' for a line with dots
% title('Light Fall-Off for Sensitivity Matrix A');
% xlabel('Source-Detector Distance (mm)');
% ylabel('Mean Signal Intensity');
% grid on;
figure;
subplot(1,2,1)
semilogy(A(sourceDetectorDistances < 30,:)')
title("Data Traces, source-detector distance < 30 mm ")
ylabel('Intensity, \muW'), xlabel('t, samples')

subplot(1,2,2)
semilogy(sourceDetectorDistances,meanIntensity,'.')
hold on
semilogy(info.pairs.r3d,mean(data,2),'.')
title("Light Falloff")
xlabel("source-detector distance, mm")
ylabel("Mean Signal Intensity")
legend("simulated","measured")
