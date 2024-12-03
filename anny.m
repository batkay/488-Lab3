%% HRF Simulation 
load("hrf_DOT3.mat");
load("NeuroDOT_Data_Sample_CCW1.mat")

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


% %% HRF simulation
% 
% % Define on/off times and repetitions
% t_on = 9;
% t_off = 36 - t_on; 
% num_rep = 4; % Number of corners
% 
% % Create timeline
% t = 1:num_rep * (t_on + t_off);
% 
% % Define HRF and convolve with stimuli
% stimuli = ones(t_on, 1); % Stimulus lasts 't_on' seconds
% data_h = conv(hrf, stimuli, 'full'); % HRF-convolved stimulus response
% 
% % Define spatial source positions in voxel coordinates
% source_positions = [
%     5, 5, 5;          % Corner 1
%     nX-5, 5, 5;       % Corner 2
%     5, nY-5, 5;       % Corner 3
%     nX-5, nY-5, 5     % Corner 4
% ];
% 
% x = zeros(nX * nY * nZ, length(t)); % Temporal data for all voxels
% y = zeros(size(voxCrd, 1), 1);      
% xVisual = zeros([nX nY nZ length(t)]);
% 
% % Convert source positions to linear indices
% source_indices = sub2ind([nX, nY, nZ], ...
%     source_positions(:,1), source_positions(:,2), source_positions(:,3));
% 
% 
% % Initialize start and end time storage
% start_t = zeros(num_rep, size(source_positions, 1));
% end_t = zeros(num_rep, size(source_positions, 1));
% 
% for rep = 1:num_rep
%     for loc = 1:size(source_positions, 1)
%         % Calculate start and end times
%         st = (rep - 1) * num_rep * t_on + (loc - 1) * (t_on) + 1;
%         et = st + length(data_h) - 1;
% 
%         % Clip end time to not exceed timeline length
%         if et > length(t)
%             et = length(t);
%         end
%         start_t(rep, loc) = st;
%         end_t(rep, loc) = et;
% 
%         % Apply HRF-convolved stimulus at the corresponding index
%         idx = source_indices(loc); % Linear index for voxel
%         x(idx, st:et) = data_h(1:(et - st + 1));
%     end
% end
% 
% 
% %% Visualization
% 
% % Loop through all time points
% for frame = 1:length(t)
%     % Extract the current frame for visualization
%     % Assuming you want to visualize the slice at z = 5
%     current_frame = reshape(x(:, frame), [nX, nY, nZ]); % Reshape to 3D volume
%     current_slice = current_frame(:, :, 5); % Extract 2D slice at z = 5
% 
%     % Visualize the current slice
%     imagesc(current_slice, [0, max(x(:))]); % Normalize the colormap
%     colorbar; % Add a colorbar for scale
%     title(['Stimulus Progression at Time = ', num2str(frame)]);
%     xlabel('X Position (mm)');
%     ylabel('Y Position (mm)');
% 
%     % Pause briefly to create an animation effect
%     pause(0.05);
% end


%% Inverse problem up on the hrf analysis
% lambda initialization:
[U, S, V] = svd(A, 'econ');
largest_singular_value = S(1,1);  % The largest singular value is the first element in S

% Set lambda as a fraction of the largest singular value
lambda = largest_singular_value / 1;  % This sets lambda to 1% of the largest singular value
A_reg = A'*(A*A'+lambda*eye(size(A,1)))^-1;
% depth = 5;% extract the slice at the DEPTH
% figure;
% 
% % Loop through all time points
% for frame = 1:length(t)
%     y = A*x(:, frame);
%     x_est = A_reg *y;
%     current_frame = reshape(x_est, [nX, nY, nZ]); 
%     current_slice = current_frame(:, :, depth); 
% 
%     % Visualize the current slice
%     imagesc(current_slice); % Normalize the colormap
%     colorbar; 
%     title(['Stimulus Progression at Time = ', num2str(frame),' at depth =',num2str(depth)]);
%     xlabel('X Position (mm)');
%     ylabel('Y Position (mm)');
% 
%     pause(0.05);
% end

%% Reconstruction of the image
% 
% data_wave_1 = data(1:round(end/2),:);
% depth = 1;
% for frame = 1:size(data_wave_1,2)
%     x_est = A_reg *data_wave_1(:, frame);
%     current_frame = reshape(x_est, [nX, nY, nZ]); 
%     current_slice = current_frame(:, :, depth); 
% 
%     % Visualize the current slice
%     imagesc(current_slice,[0, max(data_wave_1(:))]); % Normalize the colormap
%     colorbar; 
%     title(['Stimulus Progression at Time = ', num2str(frame),' at depth =',num2str(depth)]);
%     xlabel('X Position (mm)');
%     ylabel('Y Position (mm)');
% 
%     pause(0.05);
% end

% 
%% STD filtering

data_std = data(1:round(end/2), :);
std_values = std(data_std, 0, 2); % Compute standard deviation along rows
threshold = 0.00013; % Define the threshold
valid_rows = std_values <= threshold; % Identify rows with std <= threshold
data_wave_1 = data_std(valid_rows, :); % Retain only valid rows
A_reg_seg =  A_reg(:,valid_rows);

depth_range = 1:nZ; % Specify the range of depths
pause_time = 0.1; % Time to pause between frames
frame= size(data_std,2); % with range of size(data_wave_1, 2)

y = -log(bsxfun(@times,data_wave_1,1./mean(data_wave_1,2)));
x_est_t = zeros(nX*nY*nZ,frame);
% for depth = depth_range
for frame = 1:size(data_wave_1, 2)
        x_est = A_reg_seg * y(:, frame);
        x_est_t(:,frame)=x_est;
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
% synchpts = info.paradigm.synchpts(info.paradigm.Pulse_2);
% x_seg = x_est_t(:, synchpts(5):synchpts(5+1));
filtered = bandpass(x_est_t', [0.01, 0.06])';
x_est_t_reshaped = reshape(filtered, [nX, nY, nZ, size(x_est_t,2)]);

figure;
sliceViewer(squeeze(x_est_t_reshaped(:,:,10,:)),"Colormap",hot(256)); % Remove singleton dimensions
colorbar;
