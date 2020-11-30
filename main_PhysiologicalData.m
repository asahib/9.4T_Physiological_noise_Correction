close all, clear all

%Physiological data directory
%addpath('/home/jloureiro/Desktop/treated_physio');
addpath('/home/joana/Desktop/Data/fMRIStab05/fMRIStab05Physio');

%load physiological data
%matfiles = dir(fullfile('/home/jloureiro/Desktop/treated_physio'));
matfiles = dir(fullfile('/home/joana/Desktop/Data/fMRIStab05/fMRIStab05Physio'));
%addpath('/home/jloureiro/Desktop');
[DetrendedDataTotal, trigger_out_Total, PhysioData_trigger_out, rawData, mintabTotal]=DetrendedFINAL(matfiles);


% %Frequency analysis of the EPIs--------------------------------------------
% %fMRIStab05:see pg.44 - notebook2 - to see the correspondences of the volumes.
% V1=spm_read_vols(spm_vol(spm_select));%nr.12 fMRIStab05
% V2=spm_read_vols(spm_vol(spm_select));%nr.16 fMRIStab05
% V3=spm_read_vols(spm_vol(spm_select));%nr.18 fMRIStab05
% V4=spm_read_vols(spm_vol(spm_select));%nr.14 fMRIStab05 
% V5=spm_read_vols(spm_vol(spm_select));%nr.20 fMRIStab05
% VTotal={V1;V2;V3;V4;V5};
% 
% % Avgpower_psd=zeros(:,1);
% % Avgpower_psd_resp=zeros(:,1);
% % Avgpower_psd_pulse=zeros(:,1);
% Avgpower_psd_total=cell(length(VTotal));
% Avgpower_psd_resp_total=cell(length(VTotal));
% Avgpower_psd_pulse_total=cell(length(VTotal));
% 
% % DetermineVar_resp=zeros(size(Avgpower_psd,1), size(Avgpower_psd,2),2);
% % DetermineVar_pulse=zeros(size(Avgpower_psd,1), size(Avgpower_psd,2),2);
% % DetermineVar_resp_total=cell(length(VTotal));
% % DetermineVar_pulse_total=cell(length(VTotal));
% % VarResp_total=cell(length(VTotal));
% % VarPulse_total=cell(length(VTotal));
% 
% for vol=1:length(VTotal)
%     Vgeral=VTotal{vol};
%     nfft=2^nextpow2(length(Vgeral));
%     DetermineVar_resp=[];
%     DetermineVar_pulse=[];
%     for slice=1:size(Vgeral,3)
%         for i=1:size(Vgeral,1)
%             for j=1:size(Vgeral,2)
%                 Vol=squeeze(Vgeral(i,j,slice,:))';
%                 FVol=fft(Vol,nfft)/length(Vol);
%                 abs_FVol=abs(FVol(1:nfft/2+1)); 
%                 Hpsd = dspdata.psd(abs_FVol,'Fs',1000);
%                 Avgpower_psd(i,j,slice)=avgpower(Hpsd);
%                 Avgpower_psd_resp(i,j,slice)=avgpower(Hpsd, [1.1 1.2]);
%                 Avgpower_psd_pulse(i,j,slice)=avgpower(Hpsd, [0.25 0.35]);
%             end 
%         end 
%     end 
%     Avgpower_psd_total{vol}=Avgpower_psd;
%     Avgpower_psd_resp_total{vol}=Avgpower_psd_resp;
%     Avgpower_psd_pulse_total{vol}=Avgpower_psd_pulse;
%     
%     
%     %Determine the percentage of noise due to respiration and pulsefrom the
%     %fourier analysis
% 
%     DetermineVar_resp(:,:,1)=Avgpower_psd(:,:,2);
%     DetermineVar_resp(:,:,2)=Avgpower_psd_resp(:,:,2);
%     DetermineVar_pulse(:,:,1)=Avgpower_psd(:,:,2);
%     DetermineVar_pulse(:,:,2)=Avgpower_psd_pulse(:,:,2);
%     
%     DetermineVar_resp_total{vol}=DetermineVar_resp;
%     DetermineVar_pulse_total{vol}=DetermineVar_pulse;
%     var_resp=var(DetermineVar_resp,0,3);
%     var_pulse=var(DetermineVar_pulse,0,3);
%     VarResp_total{vol}=var_resp;
%     VarPulse_total{vol}=var_pulse;
% 
% end 

%% Plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Ploting raw physio with mintab vector to check

figure(1)
subplot(2,4,1)
plot(trigger_out_Total{1})
hold on
plot(mintabTotal{1,1},mintabTotal{1,2}, '--gs')
title('Pulse MZ_4sl_0.9_TE20_TR300_P4 (12,13)');

subplot(2,4,5)
plot(trigger_out_Total{2})
hold on
plot(mintabTotal{2,1},mintabTotal{2,2}, '--gs')
title('respiration MZ_4sl_0.9_TE20_TR300_P4 (12, 13)');

subplot(2,4,2)
plot(trigger_out_Total{3})
hold on
plot(mintabTotal{3,1},mintabTotal{3,2}, '--gs')
title('Pulse MZ_4sl_0.9_TE20_TR300_P4 (14,15)');

subplot(2,4,6)
plot(trigger_out_Total{4})
hold on
plot(mintabTotal{4,1},mintabTotal{4,2}, '--gs')
title('respiration MZ_4sl_0.9_TE20_TR300_P4 (14,15)');

subplot(2,4,3)
plot(trigger_out_Total{5})
hold on
plot(mintabTotal{5,1},mintabTotal{5,2}, '--gs')
title('Pulse MZ_4sl_0.9_TE20_TR1000_P4 (16,17)');

subplot(2,4,7)
plot(trigger_out_Total{6})
hold on
plot(mintabTotal{6,1},mintabTotal{6,2}, '--gs')
title('respiration MZ_4sl_0.9_TE20_TR1000_P4 (16,17)');

subplot(2,4,4)
plot(trigger_out_Total{7})
hold on
plot(mintabTotal{7,1},mintabTotal{7,2}, '--gs')
title('Pulse MZ_4sl_0.9_TE20_TR300_P4 (20,21)');

subplot(2,4,8)
plot(trigger_out_Total{8})
hold on
plot(mintabTotal{8,1},mintabTotal{8,2}, '--gs')
title('respiration MZ_4sl_0.9_TE20_TR300_P4 (20,21)');

%Ploting the detrended data to check
%t=(0:L-1)*T;

figure(2)
subplot(2,4,1)
plot(DetrendedDataTotal{1})
hold on
plot(trigger_out_Total{1}, '--g')
title('Detrended Pulse MZ_4sl_0.9_TE20_TR300_P4 (12,13)');

subplot(2,4,5)
plot(DetrendedDataTotal{2})
hold on
plot(trigger_out_Total{2}, '--g')
title('Detrended Resp MZ_4sl_0.9_TE20_TR300_P4 (12,13)');

subplot(2,4,2)
plot(DetrendedDataTotal{3})
hold on
plot(trigger_out_Total{3}, '--g')
title('Detrended Pulse MZ 4sl 0.9 TE20 TR300 P4 (nr.14,15)');

subplot(2,4,6)
plot(DetrendedDataTotal{4})
hold on
plot(trigger_out_Total{4}, '--g')
title('Detrended Resp. MZ 4sl 0.9 TE20 TR300 P4 (nr.14,15)');

subplot(2,4,3)
plot(DetrendedDataTotal{5})
hold on
plot(trigger_out_Total{5}, '--g')
title('Detrended Pulse MZ 4sl 0.9 TE20 TR1000 P4 (nr.16,17)');

subplot(2,4,7)
plot(DetrendedDataTotal{6})
hold on
plot(trigger_out_Total{6}, '--g')
title('Detrended Resp. MZ 4sl 0.9 TE20 TR1000 P4 (nr.16,17)');

subplot(2,4,4)
plot(DetrendedDataTotal{7})
hold on
plot(trigger_out_Total{7}, '--g')
title('Detrended Pulse MZ 4sl 0.9 TE20 TR300 P4 (nr.20,21)');

subplot(2,4,8)
plot(DetrendedDataTotal{8})
hold on
plot(trigger_out_Total{8}, '--g')
title('Detrended Pulse MZ 4sl 0.9 TE20 TR300 P4 (nr.20,21)');

% %Frequency Analysis of the physiological data
% 
% figure(3)
% Fs=50;
% T=1/Fs;
% 
% subplot(4,2,1)
% L=length(DetrendedDataTotal{1});
% nfft=2^nextpow2(L);
% Y=fft(DetrendedDataTotal{1},nfft)/L;
% f=Fs/2*linspace(0,1,nfft/2+1);
% plot(f, 2*abs(Y(1:nfft/2+1)));
% 
% subplot(4,2,2)
% L=length(DetrendedDataTotal{2});
% nfft=2^nextpow2(L);
% Y=fft(DetrendedDataTotal{2},nfft)/L;
% f=Fs/2*linspace(0,1,nfft/2+1);
% plot(f, 2*abs(Y(1:nfft/2+1)));
% 
% subplot(4,2,3)
% L=length(DetrendedDataTotal{3});
% nfft=2^nextpow2(L);
% Y=fft(DetrendedDataTotal{3},nfft)/L;
% f=Fs/2*linspace(0,1,nfft/2+1);
% plot(f, 2*abs(Y(1:nfft/2+1)));
% 
% subplot(4,2,4)
% L=length(DetrendedDataTotal{4});
% nfft=2^nextpow2(L);
% Y=fft(DetrendedDataTotal{4},nfft)/L;
% f=Fs/2*linspace(0,1,nfft/2+1);
% plot(f, 2*abs(Y(1:nfft/2+1)));
% 
% subplot(4,2,5)
% L=length(DetrendedDataTotal{5});
% nfft=2^nextpow2(L);
% Y=fft(DetrendedDataTotal{5},nfft)/L;
% f=Fs/2*linspace(0,1,nfft/2+1);
% plot(f, 2*abs(Y(1:nfft/2+1)));
% 
% subplot(4,2,6)
% L=length(DetrendedDataTotal{6});
% nfft=2^nextpow2(L);
% Y=fft(DetrendedDataTotal{6},nfft)/L;
% f=Fs/2*linspace(0,1,nfft/2+1);
% plot(f, 2*abs(Y(1:nfft/2+1)));
% 
% subplot(4,2,7)
% L=length(DetrendedDataTotal{7});
% nfft=2^nextpow2(L);
% Y=fft(DetrendedDataTotal{7},nfft)/L;
% f=Fs/2*linspace(0,1,nfft/2+1);
% plot(f, 2*abs(Y(1:nfft/2+1)));
% 
% subplot(4,2,8)
% L=length(DetrendedDataTotal{8});
% nfft=2^nextpow2(L);
% Y=fft(DetrendedDataTotal{8},nfft)/L;
% f=Fs/2*linspace(0,1,nfft/2+1);
% plot(f, 2*abs(Y(1:nfft/2+1)));
% 
% %Plots for the Frequency analysis of the EPIs------------------------------
% %transversal sequences-------------------------------
% figure(4)
% 
% subplot(3,3,1)
% Avg1=Avgpower_psd_total{1};
% Avg_resp1=Avgpower_psd_resp_total{1};
% Avg_pulse1=Avgpower_psd_pulse_total{1};
% imagesc(rot90(Avg1(:,:,2)));
% %title('Whole Spectrum of frequencies')
% subplot(3,3,4)
% imagesc(rot90(Avg_resp1(:,:,2)));
% %title('Respiration frequency range')
%  subplot(3,3,7)
% imagesc(rot90(Avg_pulse1(:,:,2)));
% %title('Pulse frequency range')
% 
% subplot(3,3,3)
% Avg2=Avgpower_psd_total{2};
% Avg_resp2=Avgpower_psd_resp_total{2};
% Avg_pulse2=Avgpower_psd_pulse_total{2};
% imagesc(rot90(Avg2(:,:,2)));
% %title('Whole Spectrum of frequencies')
% subplot(3,3,6)
% imagesc(rot90(Avg_resp2(:,:,2)));
% %title('Respiration frequency range')
%  subplot(3,3,9)
% imagesc(rot90(Avg_pulse2(:,:,2)));
% %title('Pulse frequency range')
% 
% subplot(3,3,2)
% Avg3=Avgpower_psd_total{3};
% Avg_resp3=Avgpower_psd_resp_total{3};
% Avg_pulse3=Avgpower_psd_pulse_total{3};
% imagesc(rot90(Avg3(:,:,2)));
% %title('Whole Spectrum of frequencies')
% subplot(3,3,5)
% imagesc(rot90(Avg_resp3(:,:,2)));
% %title('Respiration frequency range')
%  subplot(3,3,8)
% imagesc(rot90(Avg_pulse3(:,:,2)));
% %title('Pulse frequency range')
% 
% %coronal sequences------------------------------------------
% 
% figure(5)
% 
% subplot(3,2,1)
% Avg4=Avgpower_psd_total{4};
% Avg_resp4=Avgpower_psd_resp_total{4};
% Avg_pulse4=Avgpower_psd_pulse_total{4};
% imagesc(rot90(Avg4(:,:,2)));
% title('Whole Spectrum of frequencies')
% subplot(3,2,3)
% imagesc(rot90(Avg_resp4(:,:,2)));
% title('Respiration frequency range')
%  subplot(3,2,5)
% imagesc(rot90(Avg_pulse4(:,:,2)));
% title('Pulse frequency range')
% 
% subplot(3,2,2)
% Avg5=Avgpower_psd_total{5};
% Avg_resp5=Avgpower_psd_resp_total{5};
% Avg_pulse5=Avgpower_psd_pulse_total{5};
% imagesc(rot90(Avg5(:,:,2)));
% title('Whole Spectrum of frequencies')
% subplot(3,2,4)
% imagesc(rot90(Avg_resp5(:,:,2)));
% title('Respiration frequency range')
%  subplot(3,2,6)
% imagesc(rot90(Avg_pulse5(:,:,2)));
% title('Pulse frequency range')
% 
% %Plots for the Analysis of variance EPIs-----------------------------------
% 
% %transversal sequences--------------
% 
% figure(6)
% 
% subplot(2,3,1)
% VarResp1=VarResp_total{1};
% VarPulse1=VarPulse_total{1};
% imagesc(rot90(VarResp1(:,:)));
% title('Whole Spectrum of frequencies')
% subplot(2,3,4)
% imagesc(rot90(VarPulse1(:,:)));
% title('Respiration frequency range')
%  
% subplot(2,3,2)
% VarResp2=VarResp_total{2};
% VarPulse2=VarPulse_total{2};
% imagesc(rot90(VarResp2(:,:)));
% title('Whole Spectrum of frequencies')
% subplot(2,3,5)
% imagesc(rot90(VarPulse2(:,:)));
% title('Respiration frequency range')
% 
% subplot(2,3,3)
% VarResp3=VarResp_total{3};
% VarPulse3=VarPulse_total{3};
% imagesc(rot90(VarResp3(:,:)));
% title('Whole Spectrum of frequencies')
% subplot(2,3,6)
% imagesc(rot90(VarPulse3(:,:)));
% title('Respiration frequency range')
% 
% 
% %coronal sequences-----------------------
% 
% figure(7)
% 
% subplot(2,2,1)
% VarResp4=VarResp_total{4};
% VarPulse4=VarPulse_total{4};
% imagesc(rot90(VarResp1(:,:)));
% title('Whole Spectrum of frequencies')
% subplot(2,3,3)
% imagesc(rot90(VarPulse4(:,:)));
% title('Respiration frequency range')
%  
% subplot(2,2,2)
% VarResp5=VarResp_total{5};
% VarPulse5=VarPulse_total{5};
% imagesc(rot90(VarResp5(:,:)));
% title('Whole Spectrum of frequencies')
% subplot(2,2,4)
% imagesc(rot90(VarPulse5(:,:)));
% title('Respiration frequency range')
% 


% subplot(3,2,2)
% imagesc(rot90(Avgpower_psd_noNoise(:,:,2)));
% title('Without including physiological noise frequency range')
% subplot(3,2,4)
% imagesc(rot90(Avgpower_psd_noResp(:,:,2)));
% title('Without respiration frequency range ')
% subplot(3,2,6)
% imagesc(rot90(Avgpower_psd(:,:,2)));
% title('Without pulse frequency range')



