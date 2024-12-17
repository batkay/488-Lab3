%% Inverse_problem

% getting x from Ay

%% Diffuse optics Green's function for an infinite medium, computed in a rectangular slab
load('Data/NeuroDOT_Data_Sample_OUT1.mat'); % data, info, flags


mua = .0192; % flags.op.mua_gray=[0.0180,0.0192];
musp = 0.6726; % flags.op.musp_gray=[0.8359,0.6726];
nu = 1.4; 
D = 1/(3*(mua+musp)); % diffusion coefficient
mu_eff = sqrt(mua/D); % effective attenuation coefficien

% Define bounds on medium
xBnds = [-75 75]; yBnds = [-45 45]; zBnds = [1 30];  
mmX = 2; mmY = 2; mmZ = 2; 

[Y X Z] = meshgrid(yBnds(1):mmY:yBnds(2), xBnds(1):mmX:xBnds(2), zBnds(1):mmZ:zBnds(2)); % generate coordinates for slab 

voxCrd = double([X(:) Y(:) Z(:)]); % coordinates (mm) for each voxel, reshaped as 1D vector

nX = size(X,1); nY = size(X,2); nZ = size(X,3); % volume is nX x nY x nZ voxels
srcPos =[info.optodes.spos2, ones(size(info.optodes.spos2, 1), 1)];
detPos = [info.optodes.dpos2, ones(size(info.optodes.dpos2, 1), 1)];


%% A Matrix calculation
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
        % Gs = 1./(4*pi*D*r2).*exp(-mu_eff*r2);
        % ./Gs;
         A(itr,:) = A_tmp;
        itr =itr+1;
    end
end

Gs = reshape(Gs, size(Gs, 1) * size(Gs, 2), []);
A = A ./ Gs;

[U, S, V] = svd(A, 'econ');
largest_singular_value = S(1,1);  % The largest singular value is the first element in S

% Set lambda as a fraction of the largest singular value
lambda = largest_singular_value / 10;  % This sets lambda to 1% of the largest singular value


% lambda = .01;
% lambda = 0.01;


%%
% frameSize = info.io.nframe;

% y = -log(data./mean(data')');
y = log(bsxfun(@times,data,1./mean(data,2)));

% out = y(1:size(y, 1)/2, :);
out = y;
out = out( (info.pairs.NN==2 | info.pairs.NN == 1) & info.pairs.WL == 2, :);
A_seg = A((info.pairs.NN(info.pairs.WL == 2) == 2 | info.pairs.NN(info.pairs.WL == 2) == 1), :);

deviation = std(out, 0, 2);
threshold = 0.1;
valid_rows = deviation <= threshold; % Identify rows with std <= threshold
data_wave_1 = out(valid_rows, :); % Retain only valid rows
A_seg =  A_seg(valid_rows, :);

valid_rows = max(abs(data_wave_1), [], 2) <= threshold; % Identify rows with std <= threshold
data_wave_1 = data_wave_1(valid_rows, :); % Retain only valid rows
A_seg =  A_seg(valid_rows, :);

% for i = 1:size(out, 1)
%     % if(abs(max(filtered(i, :))) > 0.2)
%     %     filtered(i, :) = 0;
%     % end
% 
%     if (deviation(i) > 0.13)
%         out(i, :) = 0;
%     end
% end


% [filtered, D] = bandpass(data_wave_1', [0.025, .5], 11, StopbandAttenuation=60);
% 
% filtered = filtered';


[B, BatA] = butter(2, [.02, 1]/(11/2), "bandpass"); %0.02, 0.03
filtered = filter(B, BatA, data_wave_1);
filtered = movmean(filtered', 10)';

averageDat = zeros(size(filtered, 1), 400, 11);

for i = 1:length(info.paradigm.Pulse_2)
    idx1 = info.paradigm.synchpts(info.paradigm.Pulse_2(i));
    % idx2 = info.paradigm.synchpts(info.paradigm.Pulse_2(i)) + 400;

    averageDat(:, :, i) = [filtered(:, idx1:idx1+399)];
end

average = mean(averageDat, 3);
filtered = average;
% filtered = bandpass(filtered', [1, 5], 11)';


% filtered(abs(filtered) < 0.002) = 0.00;
% [filtered, D] = bandpass(out', [1, 5], 11);

% filtered = highpass(filtered', 3, 11)';
% filtered = highpass(y', 0.1, 11)';
% filtered = bandpass(y', [0.05, 0.2], 11)';
% filtered(filtered > 0.02) = 0.2;
%filtered(abs(filtered) < 0.002) = 0.00;
% y = -log(bsxfun(@times,filtered,1./mean(filtered,2)));
%filtered = data;
% filtered = bandpass(y', [0.1, 5], 11)';
% filtered = lowpass(y', 0.1, 11)';

% y = -ln(data/mean(data')');

% y is percent change, can cutoff temporal variance
% want to remove bad y measurements
% can kill low level signals less than like 25% of max

% lowpass for 1 sorta does stuff

% samples by source detector pairs
A_reg = A_seg'*(A_seg*A_seg'+lambda*eye(size(A_seg,1)))^-1;

xapprox = A_reg * filtered;
time = size(xapprox, 2);
% M(time) = struct('cdata',[],'colormap',[]);
movie = zeros(nX, nY, time);
for t = 1:time
    % xapprox = A_reg * filtered(1:size(filtered, 1)/2, t);
    % frame = reshape(xapprox(:, 1), nX, nY, nZ);
    frame = reshape(xapprox(:, t),nX,nY,nZ);
    movie(:, :, t) = frame(:, :, 5);
end

figure;
sliceViewer(movie, "Colormap", hot(256));

% figure;
% for t = 1:time/11
%     % frame = reshape(xapprox(:, t),nX,nY,nZ);
%     % f = im2frame(frame(:, :, 5));
%     % f = im2frame(imagesc(frame(:, :, 5)));
%     % imagesc(frame(:, :, 5));
%     imagesc(movie(:, :, t), [-.003, .003]);
%     drawnow;
%     pause(0.05);
% 
%     % f = getframe();
%     % M(t) = f; 
% end

% movie(M, 1, 10);
