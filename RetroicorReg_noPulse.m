%% RETROICOR method with the physio data treated my way

function rRetroicor=RetroicorReg_noPulse(resp)

step=length(resp)/600;


TR_start=round(1:step:length(resp));
TR_start=TR_start';
numTR =length(TR_start);
TRlen = mean(diff(TR_start));
dt = 1/49.8;
sr = 1/49.8;


% % compute cardiac phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% ns=length(Physio);
% smp = 1:ns;
% 
% max_pulse=local_max(pulse);
% p_peaks=smp(max_pulse);
% pulse_mean=mean(pulse);
% 
% 
% rmpk1 = find(pulse(p_peaks) < pulse_mean);
% rmpk = rmpk1;
% p_peaks = setdiff(p_peaks,p_peaks(rmpk));
% 
% for i=1:length(pulse)
%     if i < p_peaks(1)
%         c_phs(i) = NaN;
%     elseif i >= p_peaks(end)
%         c_phs(i) = NaN;
%     else
%         prev_peak = max(find(p_peaks <=i));
%         t1 = p_peaks(prev_peak);
%         t2 = p_peaks(prev_peak+1);
%         c_phs(i) = 2*pi*(i - t1)/(t2-t1);
%     end
% end
% 
% figure
% plot(c_phs,'m');
% xlabel('Samples');
% title('Computing Cardiac Phase');
% ylabel('mV');
% 
% %indexes where there each timecourse acquisition starts - TR=1s so the
% %acquisition happens every 49.8 points of the physio data because the
% %fs(pulse)=49.8Hz
% 
% 
% TR_phs(:,1) = c_phs(TR_start);
% 
% for i = 1:2
%     dm_c_phs(:,(i*2)-1) = cos(i*TR_phs(:,1));
%     dm_c_phs(:,i*2) = sin(i*TR_phs(:,1));
% end
% 
% figure;
% imagesc(dm_c_phs);
% colormap(gray);
    
% compute Respiratory phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%normalize to range of 0 to 1
resp_range = range(resp(:));
resp_norm = zeros(length(resp));
resp_norm = (resp(:) - min(resp(:))) / resp_range;

figure;
subplot(3,1,1);
%plot(data(:,resp_col));
plot(resp_norm);
hold; 
  
%histogram-equalized transfer function between respiratory amplitude and
%resp phase

nbins = 100;
[resp_hist,bins] = hist(resp_norm,nbins);
resp_transfer_func = [0 (cumsum(resp_hist) / sum(resp_hist))];
kern_size = round(1/dt - 1);
resp_smooth = conv(resp_norm,ones(kern_size,1),'same'); %smoothed version for taking derivative
resp_diff = [diff(resp_smooth);0]; %derivative dR/dt
r_phs = pi*resp_transfer_func(round(resp_norm * nbins)+1)' .* sign(resp_diff);

plot(resp_smooth / max(resp_smooth),'g'); %plot smoothed version
axis([0 length(resp_norm) 0 1]);

subplot(3,1,2);
plot(r_phs,'m');
hold;
 for i=1:numTR
    plot([TR_start(i) TR_start(i)],[-pi pi],'g');
 end
axis([0 length(r_phs) -pi pi]);
%get TR phase
TR_phs(:,2) = r_phs(TR_start);
subplot(3,1,3);
scatter(TR_phs(:,2),resp_norm(TR_start));
axis([-pi pi 0 1]);
xlabel('Phase in Respiratory Cycle (Radians)');
ylabel('Normalized Respiration Belt (norm V) ');

    %fit Xth order fourier series to estimate phase
    order=2;
for i = 1:order
    dm_r_phs(:,(i*2)-1) = cos(i*TR_phs(:,2));
    dm_r_phs(:,i*2) = sin(i*TR_phs(:,2));
end

% Assign output matrix and Regress it out from the image

rRetroicor = [dm_r_phs];
end


