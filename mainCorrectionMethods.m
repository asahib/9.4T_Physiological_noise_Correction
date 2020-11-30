%% main for the Analysis of the diferent noise correction methods

MAG=spm_read_vols(spm_vol(spm_select));
PH=spm_read_vols(spm_vol(spm_select));

V08BOLD04=spm_read_vols(spm_vol(spm_select));
P09BOLD04=spm_read_vols(spm_vol(spm_select));
PH09BOLD04=spm_read_vols(spm_vol(spm_select));

V08BOLD05=spm_read_vols(spm_vol(spm_select));
V09BOLD05=spm_read_vols(spm_vol(spm_select));
PH09BOLD05=spm_read_vols(spm_vol(spm_select));
V1BOLD05=spm_read_vols(spm_vol(spm_select));

V08BOLD06=spm_read_vols(spm_vol(spm_select));
V09BOLD06=spm_read_vols(spm_vol(spm_select));
V1BOLD06=spm_read_vols(spm_vol(spm_select));

V09_1_BOLD08=spm_read_vols(spm_vol(spm_select));
V09_2_BOLD08=spm_read_vols(spm_vol(spm_select));

V09Stab05_300=spm_read_vols(spm_vol(spm_select));
PH09Stab05_300=spm_read_vols(spm_vol(spm_select));
V09Stab05_1000=spm_read_vols(spm_vol(spm_select));
PH09Stab05_1000=spm_read_vols(spm_vol(spm_select));

MAG=V09BOLD05;
PH=PH09BOLD05;

%% Unwraping and rescaling phase images
%%2. Unwrap phase images 
[uPH, PH]=unwrapTime(PH);

%% Detrending magnitude images 

dMAG=zeros(size(MAG,1),size(MAG,2),size(MAG,3),size(MAG,4) );
detrended = detrendata(MAG);

for slice=1:size(MAG,3)
    for i=1:size(MAG,1)
        for j=1:size(MAG,2)
            dMAG(i,j,slice,:)=(squeeze(mean(MAG(i,j,slice,:)))+squeeze(detrended(i,j,slice,:)));
        end  
    end  
end

%convert new .mat corrected images in .nii images
dMAGnii.vol=dMAG(:,:,:,:);
err=MRIwrite(dMAGnii, '/home/joana/Desktop/Data/BOLD/FMRIBOLD05/Corrected_images/dMAGnii.nii', 'double');

%% Detrending Phase images 

duPH=zeros(size(PH,1),size(PH,2),size(PH,3),size(PH,4) );
detrended = detrendata(uPH);

for slice=1:size(PH,3)
    for i=1:size(PH,1)
        for j=1:size(PH,2)
            duPH(i,j,slice,:)=(squeeze(mean(uPH(i,j,slice,:)))+squeeze(detrended(i,j,slice,:)));
        end  
    end  
end

%convert new .mat corrected images in .nii images
duPHnii.vol=duPH(:,:,:,:);
err=MRIwrite(duPHnii, '/home/joana/Desktop/Data/BOLD/FMRIBOLD05/Corrected_images/duPHnii.nii', 'double');


%% Check Magnitude and the detrended Magnitude and Phase images

figure
subplot(421)
imagesc(rot90(mean(MAG(:,:,7,:),4)))
axis off
title('Raw Magnitude image');

subplot(422)
imagesc(rot90(mean(dMAG(:,:,7,:),4)));
axis off
title('Detrended Magnitude image');

subplot(423)
imagesc(rot90(mean(MAG(:,:,7,:),4)./std(MAG(:,:,7,:),[],4)));
axis off
title('tSNR map of the raw magnitude image')

subplot(424)
imagesc(rot90(mean(dMAG(:,:,7,:),4)./std(dMAG(:,:,7,:),[],4)));
axis off
title('tSNR map of the detrended magnitude image')

subplot(425)
imagesc(rot90(mean(uPH(:,:,7,:),4)))
axis off
title('Raw Magnitude image');

subplot(426)
imagesc(rot90(mean(duPH(:,:,7,:),4)));
axis off
title('Detrended Magnitude image');

subplot(427)
imagesc(rot90(mean(uPH(:,:,7,:),4)./std(uPH(:,:,7,:),[],4)));
axis off
title('tSNR map of the raw magnitude image')

subplot(428)
imagesc(rot90(mean(duPH(:,:,7,:),4)./std(duPH(:,:,7,:),[],4)));
axis off
title('tSNR map of the detrended magnitude image')

%% NewCoKMethod Implementation

%Detrending the regressor
[rNewCoK,rNewCoKdet, newMAGNewCoK]=NewCoKMethod(PH,uPH, MAG,dMAG);

%not detrending the regressor
%[rNewCoK, newMAGNewCoK]=NewCoKMethod(PH,uPH, MAG,dMAG);

%Correcting the phase images
[rNewCoKph,rNewCoKdetph, newPHNewCoK]=NewCoKMethod(PH,uPH,uPH,duPH);

%convert new MAG .mat corrected images in .nii images
newMAGNNewCoKnii.vol=newMAGNewCoK(:,:,:,:);
err=MRIwrite(newMAGNNewCoKnii, '/home/joana/Desktop/Data/BOLD/FMRIBOLD05/Corrected_images/newMAGNewCoKnii.nii', 'double');

%convert new PH .mat corrected images in .nii images
newPHNewCoKnii.vol=newMAGNewCoK(:,:,:,:);
err=MRIwrite(newPHNewCoKnii, '/home/joana/Desktop/Data/BOLD/FMRIBOLD05/Corrected_images/newPHNewCoKnii.nii', 'double');

%Fourier Analysis of a group of voxels located in the SC

%matrix that contains the rows and columns of the voxels located on the SC
%area. The first 3 voxels are on the left (of the image) SC and the last 3 on the right SC
SCvoxels=[113,105;114,104;120,106;117,114;112,115;119,113]; 
tsMAG=zeros(length(SCvoxels),size(MAG,4));
newtsMAG=zeros(length(SCvoxels),size(MAG,4));
tsPH=zeros(length(SCvoxels),size(MAG,4));
newtsPH=zeros(length(SCvoxels),size(MAG,4));

for vox=1:length(SCvoxels)
    tsMAG(vox,:)=(squeeze(dMAG(SCvoxels(vox,1),SCvoxels(vox,2),7,:)));
    tsPH(vox,:)=(squeeze(duPH(SCvoxels(vox,1),SCvoxels(vox,2),7,:)));
    betaMAG=pinvX\squeeze(tsMAG(vox,:))';
    betaPH=pinvX\squeeze(tsPH(vox,:))';
    newtsMAG(vox,:)=squeeze(tsMAG(vox,:))-betaMAG(1)*X(1,:);
    newtsPH(vox,:)=squeeze(tsPH(vox,:))-betaPH(1)*X(1,:);
end  

%plots of the Fourier analysis

figure
subplot(431)
plot(squeeze(tsMAG(1,:))), title('voxel (113,105)'), hold on, 
plot(squeeze(newtsMAG(1,:)),'r'), hold on,
plot(squeeze(tsPH(1,:)),'c'), hold on
plot(squeeze(newtsPH(1,:)),'m')
subplot(432)
plot(squeeze(tsMAG(2,:))), title('voxel (114,104)'), hold on, 
plot(squeeze(newtsMAG(2,:)),'r'), hold on,
plot(squeeze(tsPH(2,:)),'c'), hold on
plot(squeeze(newtsPH(2,:)),'m')
subplot(433)
plot(squeeze(tsMAG(3,:))), title('voxel (120,106)'), hold on, 
plot(squeeze(newtsMAG(3,:)),'r'), hold on,
plot(squeeze(tsPH(3,:)),'c'), hold on
plot(squeeze(newtsPH(3,:)),'m')

subplot(434)
plot(squeeze(abs(fft(detrend(tsMAG(1,1:150)))))), hold on, 
plot(squeeze(abs(fft(newtsMAG(1,1:150)))),'r'), hold on
plot(squeeze(abs(fft(detrend(tsPH(1,1:150))))),'c'), hold on, 
plot(squeeze(abs(fft(newtsPH(1,1:150)))),'m')
subplot(435)
plot(squeeze(abs(fft(tsMAG(2,1:150))))), hold on, 
plot(squeeze(abs(fft(newtsMAG(2,1:150)))),'r'), hold on
plot(squeeze(abs(fft(tsPH(2,1:150)))),'c'), hold on, 
plot(squeeze(abs(fft(newtsPH(2,1:150)))),'m')
subplot(436)
plot(squeeze(abs(fft(tsMAG(3,1:150))))), hold on, 
plot(squeeze(abs(fft(newtsMAG(3,1:150)))),'r'), hold on
plot(squeeze(abs(fft(tsPH(3,1:150)))),'c'), hold on, 
plot(squeeze(abs(fft(newtsPH(3,1:150)))),'m')

subplot(437)
plot(squeeze(tsMAG(4,:))), title('voxel (117,114)'), hold on, 
plot(squeeze(newtsMAG(4,:)),'r'), hold on,
plot(squeeze(tsPH(4,:)),'c'), hold on
plot(squeeze(newtsPH(4,:)),'m')
subplot(438)
plot(squeeze(tsMAG(5,:))), title('voxel (112,115),'), hold on, 
plot(squeeze(newtsMAG(5,:)),'r'), hold on,
plot(squeeze(tsPH(5,:)),'c'), hold on
plot(squeeze(newtsPH(5,:)),'m')
subplot(439)
plot(squeeze(tsMAG(6,:))), title('voxel (119,113)'), hold on, 
plot(squeeze(newtsMAG(6,:)),'r'), hold on,
plot(squeeze(tsPH(6,:)),'c'), hold on
plot(squeeze(newtsPH(6,:)),'m')

subplot(4,3,10)
plot(squeeze(abs(fft(tsMAG(4,1:150))))), hold on, 
plot(squeeze(abs(fft(newtsMAG(4,1:150)))),'r'), hold on
plot(squeeze(abs(fft(tsPH(4,1:150)))),'c'), hold on, 
plot(squeeze(abs(fft(newtsPH(4,1:150)))),'m')
subplot(4,3,11)
plot(squeeze(abs(fft(tsMAG(5,1:150))))), hold on, 
plot(squeeze(abs(fft(newtsMAG(5,1:150)))),'r'), hold on
plot(squeeze(abs(fft(tsPH(5,1:150)))),'c'), hold on, 
plot(squeeze(abs(fft(newtsPH(5,1:150)))),'m')
subplot(4,3,12)
plot(squeeze(abs(fft(tsMAG(6,1:150))))), hold on, 
plot(squeeze(abs(fft(newtsMAG(6,1:150)))),'r'), hold on
plot(squeeze(abs(fft(tsPH(6,1:150)))),'c'), hold on, 
plot(squeeze(abs(fft(newtsPH(6,1:150)))),'m')

legend('Mag.', 'Mag. corrected NeWCoK','PH. ts','PH corrected NewCoK')

figure
subplot(221)
imagesc(rot90(mean(dMAG(:,:,7,:),4)));
%colormap Gray;
axis off
title('Detrended Magnitude image');
subplot(222)
imagesc(rot90(mean(newMAGNewCoK(:,:,7,:),4)));
axis off
title('Magnituge image after New CoK correction');

subplot(223)
imagesc(rot90(mean(dMAG(:,:,7,:),4)./std(dMAG(:,:,7,:),[],4)));
axis off
title('tSNR map of the detrended magnitude image')

subplot(224)
imagesc(rot90(mean(newMAGNewCoK(:,:,7,:),4)./std(newMAGNewCoK(:,:,7,:),[],4)));
axis off
title('tSNR map of the magnitude image corrected with New CoK method')

%% CompCorMethod Implementation

[rCompCor, newMAGCompCor, tSTD_MAG, tSTDmask_MAG]=CompCorMethod(MAG,dMAG);

%convert new .mat corrected images in .nii images
newMAGNCompCornii.vol=newMAGCompCor(:,:,:,:);
err=MRIwrite(newMAGNCompCornii, '/home/joana/Desktop/Data/BOLD/FMRIBOLD05/Corrected_images/newMAGCompCornii.nii', double);

figure
subplot(121)
imagesc(rot90(squeeze(tSTD_MAG(7,:,:))));
title('tSTD Magnitude image')
axis off
subplot(122)
imagesc(rot90(squeeze(tSTDmask_MAG(7,:,:))));
title('tSTD>max(tSTD)*0.3')
axis off

figure
subplot(221)
imagesc(rot90(mean(dMAG(:,:,7,:),4)));
axis off
title('Detrended Magnitude image');

subplot(222)
imagesc(rot90(mean(newMAGCompCor(:,:,7,:),4)));
axis off
title('Magnituge image after CompCor correction');

subplot(223)
imagesc(rot90(mean(dMAG(:,:,7,:),4)./std(dMAG(:,:,7,:),[],4)));
axis off
title('tSNR map of the detrended magnitude image')

subplot(224) m
imagesc(rot90(mean(newMAGCompCor(:,:,7,:),4)./std(newMAGCompCor(:,:,7,:),[],4)));
axis off
title('tSNR map of the magnitude image corrected with CompCor method')

%% HighCor method implementation

[rHighCor, CorHighCor, newMAGHighCor, CorrMAP, SelectedVoxels]=HighCorMethod(MAG,dMAG,PH,uPH, tSTDmask_MAG);

%convert new .mat corrected images in .nii images
newMAGHighCornii.vol=newMAGHighCor(:,:,:,:);
err=MRIwrite(newMAGHighCornii, '/home/joana/Desktop/Data/BOLD/FMRIBOLD05/Corrected_images/newMAGHighCornii.nii', double);

figure
subplot(221)
imagesc(rot90(mean(dMAG(:,:,7,:),4)));
axis off
title('Detrended Magnitude image');

subplot(222)
imagesc(rot90(mean(newMAGHighCor(:,:,:),3)));
axis off
title('Magnituge image after HighCor correction');

subplot(223)
imagesc(rot90(mean(dMAG(:,:,7,:),4)./std(dMAG(:,:,7,:),[],4)));
axis off
title('tSNR map of the detrended magnitude image')

subplot(224)
imagesc(rot90(mean(newMAGHighCor(:,:,:),3)./std(newMAGHighCor(:,:,:),[],3)));
axis off
title('tSNR map of the magnitude image corrected with HighCor method')


%% Retroicor Method implementation

%load the physiological data
addpath('/home/joana/Desktop/Data/BOLD/FMRIBOLD05/Physio/Physio_cut_match');

%addpath('/media/DATALinux/9T_Data/FMRIBOLD05/physio/physio_cut_match');
Phys_data = dir(fullfile('/home/joana/Desktop/Data/BOLD/FMRIBOLD05/Physio/Physio_cut_match'));
%Phys_data = dir(fullfile('/media/DATALinux/9T_Data/FMRIBOLD05/physio/physio_cut_match'));
Phys_data(1:2)=[];
detrendedData_Total=cell(length(Phys_data));
Phys_dataF=cell(length(Phys_data));
trigger_outT=cell(length(Phys_data));
mintabT=cell(length(Phys_data));
maxtabT=cell(length(Phys_data));

%detrending to find the local maxima
for i=1:length(Phys_data)
    Phys_dataF{i}=load(Phys_data(i).name);
    PhysData_trigger_out = trigger_out(cell2mat(Phys_dataF(i)));
    trigger_outT{i}=PhysData_trigger_out;
    [maxtab, mintab]=peakdet(PhysData_trigger_out,900);
    maxtabT{i}=maxtab;
    mintabT{i}=mintab;
    detrendedData=detrendPhysioData(PhysData_trigger_out,mintab);
    detrendedData_Total{i}=detrendedData;
end 

%getting the maximum and minimums manually without detrending first

 [maxtabPulse, mintabPulse]=peakdet(cell2mat(trigger_outT(1)),750);
 [maxtabResp, mintabResp]=peakdet(cell2mat(trigger_outT(2)),900);

MaximaNew{1,1}=mintabResp;
MaximaNew{1,2}=maxtabResp;
Maximanew{2,1}=mintabPulse;
MaximaNew{2,2}=maxtabPulse;

for y=1:length(mintabResp)
    if mintabResp(y,2)==0
        mintabResp(y,:)=[];
    end
end
for y=1:length(maxtabResp)
    if maxtabResp(y,2)==0
        maxtabResp(y,:)=[];
    end
end
% for y=1:length(mintabPulse)
%     if mintabPulse(y,2)==0
%         mintabPulse(y,:)=[];
%     end
% end
for y=1:length(maxtabPulse)
    if maxtabPulse(y,2)==0
        maxtabPulse(y,:)=[];
    end
end

% checking physiological data:
figure
plot(trigger_outT{1})
hold on
plot(mintabPulse(:,1),mintabPulse(:,2), '--gs')
hold on
plot(maxtabPulse(:,1),maxtabPulse(:,2), '--rs')
title('Pulse ');
% [xpulsemin,ypulsemin]=ginput;
% [xpulsemax,ypulsemax]=ginput;
% extraMintabPulse=[xpulsemin ypulsemin];
% extraMaxtabPulse=[xpulsemax ypulsemax];

figure
plot(trigger_outT{2})
hold on
plot(mintabResp(:,1),mintabResp(:,2), '--gs')
hold on
plot(maxtabResp(:,1),maxtabResp(:,2), '--rs')
title('Respiration ');
% [xrespmin,yrespmin]=ginput;
% [xrespmax,yrespmax]=ginput;
% extraMintabResp=[xrespmin yrespmin];
% extraMaxtabResp=[xrespmax yrespmax];

%including the missing data on maxtabPulse, mintabPulse, maxtabResp and
%mintabResp
ExtraMaxima{1,1}=extraMintabResp;
ExtraMaxima{1,2}=extraMaxtabResp;
ExtraMaxima{2,1}=extraMintabPulse;
ExtraMaxima{2,2}=extraMaxtabPulse;

% %correcting the respiration maxima -- NOT WORKING
% extra=1;
% for i=1:length(mintabPulse)
%     if mintabPulse(i,1)>=extraMintabPulse(extra,1)
%         C = insertrows(mintabPulse,extraMintabPulse(extra,:),i);
%         extra=extra+1;
%     else 
%         if extra==length(extraMintabPulse)
%             break
%         else
%         extra=extra+1;
%         end 
%     end  
% end  



%Ploting the detrended data to check
%t=(0:L-1)*T;
figure
subplot(2,1,1)
plot( detrendedData_Total{1})
hold on
plot(trigger_outT{1}, '--g')
title('Detrended Pulse ');

subplot(2,1,2)
plot( detrendedData_Total{2})
hold on
plot(trigger_outT{2}, '--g')
title('Detrended Resp ');

%Apply RETROICOR method
[rRetroicor newMAGRetroicor]=RetroicorMethod(MAG,dMAG, trigger_outT, MaximaNew);
[rRetroicorPH newPHRetroicor]=RetroicorMethod(PH,duPH, trigger_outT, MaximaNew);

%convert new .mat corrected images in .nii images
newMAGRetroicornii.vol=newMAGRetroicor(:,:,:,:);
err=MRIwrite(newMAGRetroicornii, '/home/joana/Desktop/Data/BOLD/FMRIBOLD05/Corrected_images/newMAGRetroicornii.nii', double);
mnewMAGRetroicor(:,:,7)=mean(dMAG(:,:,7,:),4)-mean(newMAGRetroicor(:,:,7,:),4);
figure
subplot(221)
imagesc(rot90(mean(dMAG(:,:,7,:),4)));
axis off
title('Detrended Magnitude image');

subplot(222)
imagesc(rot90(mean(newMAGRetroicor(:,:,7,:),4)));
axis off
title('Magnituge image after Retroicor correction');

subplot(223)
imagesc(rot90(mean(dMAG(:,:,7,:),4)./std(dMAG(:,:,7,:),[],4)));
axis off
title('tSNR map of the detrended magnitude image')

subplot(224)
imagesc(rot90(mean(newMAGRetroicor(:,:,7,:),4)./std(newMAGRetroicor(:,:,7,:),[],4)));
axis off
title('tSNR map of the magnitude image corrected with Retroicor method')

%% CoK method implementation 

%load physiological data
addpath('/home/joana/Desktop/Data/BOLD/FMRIBOLD05/raw');
Raw_data = dir(fullfile('/home/joana/Desktop/Data/fMRIStab05/fMRIStab05Physio'));

[rCoK, CorCoK, newMAGCoK]=CoKMethod(MAG,dMAG,Raw_data);

%convert new .mat corrected images in .nii images
newMAGCoKnii.vol=newMAGCoK(:,:,:,:);
err=MRIwrite(newMAGCoKnii, '/home/joana/Desktop/Data/BOLD/FMRIBOLD05/Corrected_images/newMAGCoKnii.nii', double);

figure
subplot(221)
imagesc(rot90(mean(dMAG(:,:,7,:),4)));
axis off
title('Detrended Magnitude image');

subplot(222)
imagesc(rot90(mean(newMAGCoK(:,:,:),3)));
axis off
title('Magnituge image after CoK correction');

subplot(223)
imagesc(rot90(mean(dMAG(:,:,7,:),4)./std(dMAG(:,:,7,:),[],4)));
axis off
title('tSNR map of the detrended magnitude image')

subplot(224)
imagesc(rot90(mean(newMAGCoK(:,:,:),3)./std(newMAGCoK(:,:,:),[],3)));
axis off
title('tSNR map of the magnitude image corrected with CoK method')

%% Comparison between all methods

%tSNR maps of the different methods
figure
subplot(431)
imagesc(rot90(mean(dMAG(:,:,7,:),4)./std(dMAG(:,:,7,:),[],4)));
axis off
title('tSNR map of the detrended magnitude image')

subplot(432)
imagesc(rot90(mean(newMAGCoK(:,:,:),3)./std(newMAGCoK(:,:,:),[],3)));
axis off
title('tSNR map of the magnitude image corrected with CoK method')

subplot(433)
imagesc(rot90(mean(newMAGNewCoK(:,:,7,:),4)./std(newMAGNewCoK(:,:,:),[],3)));
axis off
title('tSNR map of the magnitude image corrected with New CoK method')

subplot(434)
imagesc(rot90(mean(dMAG(:,:,7,:),4)./std(dMAG(:,:,7,:),[],4)));
axis off
title('tSNR map of the detrended magnitude image')

subplot(435)
imagesc(rot90(mean(newMAGCoK(:,:,:),3)./std(newMAGCoK(:,:,:),[],3)));
axis off
title('tSNR map of the magnitude image corrected with CoK method')

subplot(436)
imagesc(rot90(mean(newMAGNewCoK(:,:,7,:),4)./std(newMAGNewCoK(:,:,:),[],3)));
axis off
title('tSNR map of the magnitude image corrected with New CoK method')

subplot(437)
imagesc(rot90(mean(newMAGRetroicor(:,:,:),3)./std(newMAGRetroicor(:,:,:),[],3)));
axis off
title('tSNR map of the magnitude image corrected with Retroicor method')

subplot(438)
imagesc(rot90(mean(newMAGCompCor(:,:,:),3)./std(newMAGCompCor(:,:,:),[],3)));
axis off
title('tSNR map of the magnitude image corrected with CompCor method')

subplot(439)
imagesc(rot90(mean(newMAGHighCor(:,:,:),3)./std(newMAGHighCor(:,:,:),[],3)));
axis off
title('tSNR map of the magnitude image corrected with HighCor method')

subplot(4,3,10)
imagesc(rot90(mean(newMAGRetroicor(:,:,:),3)./std(newMAGRetroicor(:,:,:),[],3)));
axis off
title('tSNR map of the magnitude image corrected with Retroicor method')

subplot(4,3,11)
imagesc(rot90(mean(newMAGCompCor(:,:,:),3)./std(newMAGCompCor(:,:,:),[],3)));
axis off
title('tSNR map of the magnitude image corrected with CompCor method')

subplot(4,3,12)
imagesc(rot90(mean(newMAGHighCor(:,:,:),3)./std(newMAGHighCor(:,:,:),[],3)));
axis off
title('tSNR map of the magnitude image corrected with HighCor method')

%tSNR values average for SC (across subjects and averaged for all the subjects)

%Frequency analysis of the images obtained with the different methods on the
%pulse and respiration range of frequencies#

%Analysing Phase stability (SD(PH)=1/(tSNR_MAG) - theoretically)

SDPH=std(duPH(:,:,:,:),[],4);
SDPHNewCoK=std(newPHNewCoK(:,:,:,:),[],4);
SDPHCompCor=std(newPHCompCor(:,:,:,:),[],4);
SDPHHighCor=std(newPHHighCor(:,:,:,:),[],4);
SDPHRetroicor=std(newPHRetroicor(:,:,:,:),[],4);

tSNRMAG=mean(dMAG(:,:,:,:),4)./std(dMAG(:,:,:,:),[],4);
tSNRNewMAGNewCoK=mean(newMAGNewCoK(:,:,:,:),4)./std(newMAGNewCoK(:,:,:,:),[],4);
tSNRNewMAGCompCor=mean(newMAGCompCor(:,:,:,4),3)./std(newMAGCompCor(:,:,:,:),[],4);
tSNRNewMAGHighCor=mean(newMAGNewCoK(:,:,:,4),3)./std(newMAGHighCor(:,:,:,:),[],4);
tSNRNewMAGRetroicor=mean(newMAGNewCoK(:,:,:,4),3)./std(newMAGRetroicor(:,:,:,:),[],4);

invtSNRMAG=1./tSNRMAG(:,:,:);
invtSNRMAGNewCoK=1./tSNRNewMAGNewCoK(:,:,:);
invtSNRMAGCompCor=1./tSNRNewMAGCompCor(:,:,:);
invtSNRMAGHighCor=1./tSNRNewMAGHighCor(:,:,:);
invtSNRMAGRetroicor=1./tSNRNewMAGRetroicor(:,:,:);

figure
subplot(231)
imagesc(rot90(squeeze(tSNRMAG(:,:,7))))
title('dMAG tSNR fMRIBOLD05 0004')
subplot(232)
imagesc(rot90(squeeze(SDPH(:,:,7))))
title('SD of the phase image')
subplot(233)
imagesc(rot90(squeeze(invtSNRMAG(:,:,7))))
title('1/tSNR of the MAG')
subplot(234)
imagesc(rot90(squeeze(tSNRNewMAGNewCoK(:,:,7))))
subplot(235)
imagesc(rot90(squeeze(SDPHNewCoK(:,:,7))))
subplot(236)
imagesc(rot90(squeeze(invtSNRMAGNewCoK(:,:,7))))


figure
subplot(531)
imagesc(rot90(tSNRMAG(:,:,7)))
subplot(532)
imagesc(rot90(SDPH(:,:,7)))
subplot(533)
imagesc(rot90(invtSNRMAG))

subplot(534)
imagesc(rot90(tSNRNewMAGNewCoK(:,:,7)))
subplot(535)
imagesc(rot90(SDPHNewCoK(:,:,7)))
subplot(536)
imagesc(rot90(invtSNRMAGNewCoK))

subplot(537)
imagesc(rot90(tSNRNewMAGCompCor(:,:,7)))
subplot(538)
imagesc(rot90(SDPHCompCor(:,:,7)))
subplot(539)
imagesc(rot90(invtSNRMAGCompCor))

subplot(5,3,10)
imagesc(rot90(tSNRNewMAGHighCor(:,:,7)))
subplot(5,3,11)
imagesc(rot90(SDPHHighCor(:,:,7)))
subplot(5,3,12)
imagesc(rot90(invtSNRMAGHighCor))

subplot(5,3,13)
imagesc(rot90(tSNRNewMAGRetroicor(:,:,7)))
subplot(5,3,14)
imagesc(rot90(SDPHRetroicor(:,:,7)))
subplot(5,3,15)
imagesc(rot90(invtSNRMAGRetroicor))

figure
subplot(331)
imagesc(rot90(tSNRMAG(:,:,7)))
subplot(532)
imagesc(rot90(SDPH(:,:,7)))
subplot(533)
imagesc(rot90(invtSNRMAG))

subplot(334)
imagesc(rot90(tSNRNewMAGNewCoK(:,:,7)))
subplot(535)
imagesc(rot90(SDPHNewCoK(:,:,7)))
subplot(536)
imagesc(rot90(invtSNRMAGNewCoK))


subplot(3,3,7)
imagesc(rot90(tSNRNewMAGRetroicor(:,:,7)))
subplot(5,3,8)
imagesc(rot90(SDPHRetroicor(:,:,7)))
subplot(5,3,9)
imagesc(rot90(invtSNRMAGRetroicor))


%Ideas for later:
%tSD of phase images Vs. SNR0 of magnitude images for each subject (see
%gisela paper.
%Predictive Value of different noise regressors
