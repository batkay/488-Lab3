%% Inverse_problem

% getting x from Ay

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
x = ones(size(A,2),1);
y = A*x;
lambda = 2; %max(A,[],"all")/400;
A_reg = A'*(A*A'+lambda*eye(size(A,1)))^-1;
x_est = A_reg *y;

reconstrution = reshape(x_est,nX,nY,nZ);
original = reshape(x,nX,nY,nZ);
figure;
plot(x');
hold on
plot(x_est');
ylabel("intentsity changes");
xlabel("tissue structures");
legend("x","x_estimate");

figure;
imagesc(squeeze(reconstrution(:,:,1)));
colorbar;
figure;
imagesc(squeeze(original(:,:,1)));
colorbar;

figure, sliceViewer(log10(reconstrution),'Colormap',hot(256));

figure, sliceViewer(log10(original),'Colormap',hot(256));


%% PSF

freq = 1; %hz
x = zeros(size(A,2),1);
% Define point source near the center at ~10 mm depth
source_position = [round(nX/2), round(nY/2), 5];  % z index of 5 to approximately get 10 mm

% Create empty x vector and set the source position
x = zeros(nX * nY * nZ, 1);
linear_idx = sub2ind([nX, nY, nZ], source_position(1), source_position(2), source_position(3));
x(linear_idx) = 1;  % Place the point source

y = A*x;
% Compute the singular value decomposition to get the largest singular value of A
[U, S, V] = svd(A, 'econ');
largest_singular_value = S(1,1);  % The largest singular value is the first element in S

% Set lambda as a fraction of the largest singular value
lambda = largest_singular_value / 100;  % This sets lambda to 1% of the largest singular value


A_reg = A'*(A*A'+lambda*eye(size(A,1)))^-1;
x_est = A_reg *y;

reconstrution = reshape(x_est,nX,nY,nZ);
original = reshape(x,nX,nY,nZ);
figure;
plot(x');
hold on
plot(x_est');
ylabel("intentsity changes");
xlabel("tissue structures");
legend("x","x_estimate");
% 
% figure;
% imagesc(squeeze(reconstrution(:,:,1)));
% colorbar;
% figure;
% imagesc(squeeze(original(:,:,1)));
% colorbar;

figure, sliceViewer((reconstrution),'Colormap',hot(256));
colorbar;
figure, sliceViewer((original),'Colormap',hot(256));
colorbar;
%% multi point

% Assume A, nX, nY, nZ, and lambda are already defined
source_positions = [round(nX/2), round(nY/2)-18, 5; round(nX/2), round(nY/2)+18, 5];  % Starting with 10 voxels apart

% Convert positions to linear indices
idx1 = sub2ind([nX, nY, nZ], source_positions(1,1), source_positions(1,2), source_positions(1,3));
idx2 = sub2ind([nX, nY, nZ], source_positions(2,1), source_positions(2,2), source_positions(2,3));

% Initialize point sources
x = zeros(nX * nY * nZ, 1);
x([idx1, idx2]) = 1;  % Set the two point sources

% Simulate measurement
y = A * x;

% Reconstruct
A_reg = A' * (A * A' + lambda * eye(size(A,1)))^-1;
x_est = A_reg * y;
reconstruction = reshape(x_est, nX, nY, nZ);
original = reshape(x,nX,nY,nZ);
figure, sliceViewer((reconstrution),'Colormap',hot(256));

figure, sliceViewer((original),'Colormap',hot(256));


permutedData = permute(reconstruction, [1, 3, 2]);  % This moves Y to the first slicing dimension

figure;
sliceViewer(permutedData, 'Colormap', hot(256));

% Extract intensity profile along the line connecting the two sources
profile = squeeze(reconstruction(round(nX/2), :, 5));
max_intensity = max(profile);
alf_max = max_intensity / 2;
% Visualize
figure; plot(profile);
xlabel('Voxel Position');
ylabel('Intensity');
title('Profile between two point sources');
yline(half_max, '--', sprintf('Half Max at %.3f', half_max));
% Fixing title with correct formatting
title(sprintf('PSF Profile at Position %.3f, %.3f, %.3f, Lambda = %.3f', source_positions(1,1), source_positions(1,2), 5, lambda));

% Analyze the profile to find the minimum between peaks and check for distinguishability
[min_val, min_idx] = min(profile(round(nY/2):round(nY/2)+10));
% if min_val < max(profile(idx1), profile(idx2)) / 2
%     disp('Sources are distinguishable');
% else
%     disp('Sources are not distinguishable');
% end

%% 
% Define regularization parameter lambda
[U, S, V] = svd(A, 'econ');
largest_singular_value = S(1,1);
lambda_values = [largest_singular_value , largest_singular_value / 2, largest_singular_value / 4]; % Different lambda values

% Choose positions for point sources (e.g., center and edges)
positions = round([1, nX/2, nX]);  % Modify these based on your grid dimensions

% Pre-allocate space to store FWHM results
fwhm_results = zeros(length(positions), length(lambda_values));
figure;

for l = 1:length(lambda_values)
    lambda = lambda_values(l);
    
    A_reg = A' * (A * A' + lambda * eye(size(A,1)))^-1;
    
    for p = 1:length(positions)
        % Create a point source at the specified position
        x = zeros(size(A,2), 1);
        x(positions(p)) = 1;
        
        % Simulate measurement y and reconstruct PSF
        y = A * x;
        x_est = A_reg * y;
        
        % Reshape to 3D for visualization (optional)
        reconstruction = reshape(x_est, nX, nY, nZ);
        
        profile = squeeze(reconstruction(positions(p), :, 1));
        
        % Calculate FWHM
        max_intensity = max(profile);
        half_max = max_intensity / 2;
        
        indices_above_half = find(profile >= half_max);
        
        % Calculate FWHM as the width between the first and last indices above half max
        if ~isempty(indices_above_half)
            fwhm = indices_above_half(end) - indices_above_half(1);
        else
            fwhm = NaN;  % In case no points are above half max (unlikely with good data)
        end
        
        % Store the FWHM result for this position and lambda
        fwhm_results(p, l) = fwhm;
        
        % Plot the PSF profile and half max line for verification (optional)
        
        plot(profile);
        hold on;
        yline(half_max, '--', 'text','Half Max,Position %d',positions(p));
        title(sprintf('PSF Profile at Position %d, Lambda = %.3f', positions(p), lambda));
        xlabel('Position');
        ylabel('Intensity');
        
    end
end
hold off;
% Display FWHM results
disp('FWHM Results (rows: positions, columns: lambda values):');
disp(fwhm_results);
