clear all, close all

%Physiological data directory
addpath ('/home/jloureiro/Dropbox/phd_tuebingen/PROJECT/Measurements_data/19-06-2013_2volunteers/fMRIStab_02/Physiologic_data/txtFiles')
addpath ('/home/jloureiro/Dropbox/phd_tuebingen/PROJECT/Measurements_data/19-06-2013_2volunteers/fMRIStab_03/Physiologic_data/txtFiles')


%load physiological data
%fMRIStab_02-------------------------------------

%Multi-band sequence:
pulse_mb02=load('PULSE_Physio_log_20130619_103722_2404970d-1068-41bc-93fa-fb342a173ff6.txt');
resp_mb02=load('RESP_Physio_log_20130619_103722_2404970d-1068-41bc-93fa-fb342a173ff6.txt');
%MZ sequence:
pulse_mz_02=load('PULSE_PmuSignals_058.txt');
resp_mz_02=load('RESP_PmuSignals_058.txt');

%Note: the rest of the data for this subject is not valid because the
%acquisition of the sequences was aborted....


%fMRIStab_03-------------------------------------

%Multi-band sequence:
pulse1_mb03=load('PULSE_Physio_log_20130619_123037_2de33b27-e134-49a6-a711-8b3fb95243a5.txt');
resp1_mb03=load('RESP_Physio_log_20130619_123037_2de33b27-e134-49a6-a711-8b3fb95243a5.txt');
%MZ sequence TR=300ms
pulse_mz_300ms_03=load('PULSE_PmuSignals_059.txt');
resp_mz_300ms_03=load('RESP_PmuSignals_059.txt');
%MZ sequence TR=1000ms
pulse_mz_1000ms_03=load('PULSE_PmuSignals_060.txt');
resp_mz_1000ms_03=load('RESP_PmuSignals_060.txt');

%Whole brain mb epi:
% pulse2_mb03=load('PULSE_Physio_log_20130619_123744_748424d1-1d10-4b3d-8fdd-6cd62f202155.txt');
% resp2_mb03=load('RESP_Physio_log_20130619_123744_748424d1-1d10-4b3d-8fdd-6cd62f202155.txt');
% pulse3_mb03=load('PULSE_Physio_log_20130619_123857_e6382ede-3c7e-48ab-9ee0-cc32652a1446.txt');
% resp3_mb03=load('RESP_ Physio_log_20130619_123857_e6382ede-3c7e-48ab-9ee0-cc32652a1446.txt');

%V=spm_read_vols(spm_vol(spm_select));


%PLOTS--------------------------------------------------------------------

%fMRIStab_02 Physio Data__________________________________________________

figure(1)

subplot(1,2,1)
plot(pulse_mb02)
hold on
plot(resp_mb02, 'r')
title('Physio fMRIStab_02 multi-band sequence')

subplot(1,2,2)
plot(pulse_mz_02) 
hold on
plot(resp_mz_02) 
title('Physio fMRIStab_02 MZ sequence 055')

%fMRIStab_03 Physio Data_____________________________________________

figure(2)
subplot(1,3,1)
plot(pulse1_mb03)
hold on
plot(resp1_mb03, 'r')
title('Physio fMRIStab_03 multi-band sequence1')

subplot(1,3,2)
plot(pulse_mz_300ms_03) 
hold on
plot(resp_mz_300ms_03) 
title('Physio fMRIStab_03 MZ sequence 059')

subplot(1,3,3)
plot(pulse_mz_1000ms_03) 
hold on
plot(resp_mz_1000ms_03) 
title('Physio fMRIStab_03 MZ sequence 060')


