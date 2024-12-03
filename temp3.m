%% Diffuse optics Green's function for an infinite medium, computed in a rectangular slab

mua = .0192; % flags.op.mua_gray=[0.0180,0.0192];
musp = 0.6726; % flags.op.musp_gray=[0.8359,0.6726];
nu = 1.4; 
D = 1/(3*(mua+musp)); % diffusion coefficient
mu_eff = sqrt(mua/D); % effective attenuation coefficien

% Define bounds on medium
xBnds = [-65 65]; yBnds = [-45 45]; zBnds = [1 30];  
mmX = 2; mmY = 2; mmZ = 2; 

[Y X Z] = meshgrid(yBnds(1):mmY:yBnds(2), xBnds(1):mmX:xBnds(2), zBnds(1):mmZ:zBnds(2)); % generate coordinates for slab 

voxCrd = double([X(:) Y(:) Z(:)]); % coordinates (mm) for each voxel, reshaped as 1D vector

nX = size(X,1); nY = size(X,2); nZ = size(X,3); % volume is nX x nY x nZ voxels
srcPos =info.optodes.spos2;
detPos =info.optodes.dpos2;

srcPos =[info.optodes.spos2,zeros(length(srcPos),1)];
detPos =[info.optodes.dpos2,zeros(length(detPos),1)];

%% A Matrix calculation
r = pdist2(srcPos, voxCrd); % distance from source to each voxel
GsAnalytic = 1./(4*pi*D*r).*exp(-mu_eff*r); % Green's function, source to voxel

rdet = pdist2(voxCrd, detPos);
GsDet = 1./(4*pi*D*rdet).*exp(-mu_eff*rdet); % Green's function, voxel to detector

r2 = pdist2(detPos,srcPos );
Gs = 1./(4*pi*D*r2).*exp(-mu_eff*r2);

nsrc =length(srcPos);
ndet =length(detPos);
A = zeros(nsrc*ndet,length(voxCrd));
A_tmp =[];
itr = 1;
for i = 1:nsrc
    for j = 1: ndet
        A_tmp = -GsAnalytic(i,:).* GsDet(:,j)';
        r2 = pdist2(detPos(i),srcPos(j));
        % Gs = 1./(4*pi*D*r2).*exp(-mu_eff*r2);
        % ./Gs;
         A(itr,:) = A_tmp;
        itr =itr+1;
    end
end

Gs = reshape(Gs, size(Gs, 1) * size(Gs, 2), []);
A = A ./ Gs;
% lambda initialization:
[U, S, V] = svd(A, 'econ');
largest_singular_value = S(1,1);  % The largest singular value is the first element in S

% Set lambda as a fraction of the largest singular value
lambda = largest_singular_value / 10;  % This sets lambda to 1% of the largest singular value

A_reg = A'*(A*A'+lambda*eye(size(A,1)))^-1;
depth = 5;% extract the slice at the DEPTH
%% STD filtering

data_std = data(1:round(end/2), :);
std_values = std(data_std, 0, 2); % Compute standard deviation along rows
threshold = 0.0013; % Define the threshold
valid_rows = std_values <= threshold; % Identify rows with std <= threshold
data_wave_1 = data_std(valid_rows, :); % Retain only valid rows
A_reg =  A_reg(:,valid_rows);

depth_range = 1:nZ; % Specify the range of depths
pause_time = 0.1; % Time to pause between frames
frame= 4857; % with range of size(data_wave_1, 2)
depth =(2);

y = -log(bsxfun(@times,data_wave_1,1./mean(data_wave_1,2)));
y_est = zeros(28*24,size(data_wave_1, 2));
x_est_t = zeros(nX*nY*nZ,frame);
% for depth = depth_range
for frame = 1:size(data_wave_1, 2)
        x_est = A_reg * y(:, frame);
        x_est_t(:,frame)=x_est;
        y_est(:,frame) = A * x_est;
        % current_frame = reshape(x_est, [nX, nY, nZ]);
        % current_slice = current_frame(:, :, depth);

    %     % Visualize the current slice
    %     imagesc(current_slice, [0, max(data_wave_1(:))]); % Normalize the colormap
    %     colorbar;
    %     title(['Stimulus Progression at Time = ', num2str(frame), ...
    %            ', Depth = ', num2str(depth)]);
    %     xlabel('X Position (mm)');
    %     ylabel('Y Position (mm)');
    % 
    %     pause(pause_time); % Pause for a short time to allow visualization
    % % end
end

%%


% data_std = data(1:round(end/2), :);
% std_values = std(data_std, 0, 2); % Compute standard deviation along rows
% threshold = 0.013; % Define the threshold
% valid_rows = std_values <= threshold; % Identify rows with std <= threshold
% data_wave_1 = data_std(valid_rows, :); % Retain only valid rows

% figure
% plot(fft(x_est_t)');
% filtered = bandpass(x_est_t', [0.001, 0.1])';
filtered = bandpass(x_est_t', [0.0001, 0.1])';
x_est_t_reshaped = reshape(filtered, [nX, nY, nZ, frame]);
figure;
sliceViewer(squeeze(x_est_t_reshaped(:,:,5,:)),"Colormap",hot(256)); % Remove singleton dimensions

colorbar

