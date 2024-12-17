% load("NeuroDOT_Data_Sample_CCW1.mat")
% 
% y = -log(bsxfun(@times,data,1./mean(data,2)));
% 
% filtered = highpass(y', 1, 11)';
% filtered(abs(filtered) > 0.02) = sign(filtered(abs(filtered) > 0.02)) * 0.02;
% filtered(abs(filtered) < 0.002) = 0.00;

% figure;
% subplot(2, 2, 1);
% imagesc(y), caxis([-.2, .2]); % need to axis limit
figure;
title("y");
subplot(2, 2, 2);
imagesc(filtered), caxis([-.1, .1]);
title("Filter")
subplot(2, 2, 3);
plot(1:length(y), y);
title("Original")
subplot(2, 2, 4);
plot(1:length(filtered), filtered);
title("Filter")

%% fft

% y = -log(bsxfun(@times,data,1./mean(data,2)));
% 
% % filtered = bandpass(y', [0.05, 0.2], 11)';
% filtered = y;
% for i = 1:size(filtered, 1)
%     if(abs(max(filtered(i, :))) > 0.2)
%         filtered(i, :) = 0;
%     end
% end
% filtered(abs(filtered) < 0.002) = 0.00;

fs = 11;
ff = [0:fs/size(filtered,2):fs/2]; % frequency vector 
Yf = fft(filtered,[],2); % fft along 2nd dimension of y, all measurements 


fy = [0:fs/size(y, 2):fs/2];
Y = fft(y, [], 2);

figure;
subplot(1, 2, 1);
plot(ff,mean(abs(Yf(:,1:length(ff))),1));
subplot(1, 2, 2);
plot(1:length(filtered), filtered(10, :));

figure;
subplot(1, 2, 1);
plot(fy,mean(abs(Y(:,1:length(fy))),1));
subplot(1, 2, 2);
plot(1:length(y), y(10, :));