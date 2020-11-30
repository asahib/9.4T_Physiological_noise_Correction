%% get central of k-space

[CoK, newCoK, image_obj] =GetCoK_modified_Gisela('meas_MID1078_mzBOLD_1190_4sl_0_9mSineNoDiCoTE20P4TR300_FID1703.dat');
dummy = squeeze(image_obj{210,:,54,1,:});

%% Regress the signal of one chanel out

CoKStab050012_ch20=zeros(size(VStab050012,1),size(VStab050012,2));
 newVStab050012CoKch20Corrected=zeros(size(VStab050012,1),size(VStab050012,2), 1000);
for i=1:size(VStab050012,1)
    for j=1:size(VStab050012,2)
        CoKStab050012_ch20(i,j)=regress(squeeze(VStab050012(i,j,3,1:1000)),squeeze(newCoK(20,:,3))');
        newVStab050012CoKch20Corrected(i,j,:)=mean(squeeze(VStab050012(i,j,3,1:1000)))+(squeeze(VStab050012(i,j,3,1:1000))-(CoKStab050012_ch20(i,j)*squeeze(newCoK(20,:,3)))');

    end 
end

newSDretroicor0012=std(Vretroicor0012(:,:,3,:),[],4);
tSNRretroicor=(mean(Vretroicor0012(:,:,3,:),4))./newSDretroicor0012(:,:);
%PLOTS
figure
subplot(231)
imagesc(rot90(tSNRnocor));
title('tSNR maps with no correction')
axis off
subplot(232)
imagesc(rot90(tSNR));
title('tSNR maps after CoK correction ch20')
axis off
subplot(233)
imagesc(rot90(tSNRretroicor));
title('tSNR maps after RETROICOR correction')
axis off

subplot(234)
imagesc(rot90(tSNRnocor));
axis off
subplot(235)
imagesc(rot90(tSNR));
axis off
subplot(236)
imagesc(rot90(tSNRretroicor));
axis off






figure
imagesc(rot90(tSNR))
roi=imfreehand();
binaryImage=roi.createMask();
imshow(binaryImage);
masked_tSNR=tSNR;
masked_tSNR(~binaryImage)=0;
imshow(masked_tSNR);
meantSNR=mean(masked_tSNR(binaryImage));

%% Regress CoK using SVD
nrChan=size(newCoK,1);
nrSlices=size(newCoK,3);
nrFrames=size(newCoK,2);
%t=0:0.3:300-0.3;
%n
SVDMatrixTotal=cell(nrFrames);
UTotal=cell(nrFrames);
VTotal=cell(nrFrames);
AfterSVDTotal=cell(nrFrames);
lambdaTotal=cell(nrFrames);
lambda=zeros(31, nrFrames);
firstComp=zeros(nrFrames,size(newCoK,2));
Ltime=nrFrames;

for sl=1:nrSlices
    [U,SVDMatrix,V]=svd(squeeze(newCoK(:,:,sl)));
    SVDMatrixTotal{sl}=SVDMatrix;
    UTotal{sl}=U;
    VTotal{sl}=V;
    lambda(:,sl)=diag(SVDMatrix);
    lambdaTotal{sl}=lambda;
    AfterSVD=U*SVDMatrix*V';
    AfterSVDTotal{sl}=AfterSVD;
    firstComp(sl,:)=AfterSVD(1,:);
end  
[U,SVDMatrix,V]=svd(timeSeries(:,:));
lambda=diag(SVDMatrix);
variances=lambda.*lambda;
bar(variances(1:30));
PCs=3;
VV=V(:,1:PCs);
Y=VV'*timeSeries';
XX=VV*Y;


CoKStab050012_SVDcorrected=zeros(size(VStab050012,1),size(VStab050012,2));
 newVStab050012CoKSVDCorrected=zeros(size(VStab050012,1),size(VStab050012,2), 1000);
for i=1:size(VStab050012,1)
    for j=1:size(VStab050012,2)
        CoKStab050012_SVDcorrected(i,j)=regress(squeeze(VStab050012(i,j,3,1:1000)),squeeze(firstComp(3,:))');
        newVStab050012CoKSVDCorrected(i,j,:)=mean(squeeze(VStab050012(i,j,3,1:1000)))+(squeeze(VStab050012(i,j,3,1:1000))-(CoKStab050012_ch20(i,j)*squeeze(firstComp(3,:)))');

    end 
end


%plot(100*lambda(1,:)./sum(lambda,1));%percent variance represented by the
%first eigenvector



figure(1)
for ch=1:6
    subplot(3,2,ch)
    plot(t,detrend(squeeze(newCoK(ch,:,2)))/max(detrend(squeeze(newCoK(ch,:,2)))))%slice nr.2
    title(['ch no ' num2str(ch)]);
%     %hold on
%     plot(tPulse, detrendedData/max(detrendedData),'g');
%     hold on,
%     plot(tResp,DetrendedDataTotal{2}/max(DetrendedDataTotal{2}),'r');
end

figure(2)
for ch=7:12
    subplot(3,2,ch-6)
    plot(t,detrend(squeeze(newCoK(ch,:,2)))/max(detrend(squeeze(newCoK(ch,:,2)))))%slice nr.2
    title(['ch no ' num2str(ch)]);
%     hold on
%     plot(tPulse, detrendedData/max(detrendedData),'g');
%     hold on,
%     plot(tResp,DetrendedDataTotal{2}/max(DetrendedDataTotal{2}),'r');
end

figure(3)
for ch=13:18
    subplot(3,2,ch-12)
    plot(t,detrend(squeeze(newCoK(ch,:,2)))/max(detrend(squeeze(newCoK(ch,:,2)))))%slice nr.2
    title(['ch no ' num2str(ch)]);
%     hold on
%     plot(tPulse, detrendedData/max(detrendedData),'g');
%     hold on,
%     plot(tResp,DetrendedDataTotal{2}/max(DetrendedDataTotal{2}),'r');
end

figure(4)
for ch=18:23
    subplot(3,2,ch-17)
    plot(t,detrend(squeeze(newCoK(ch,:,2)))/max(detrend(squeeze(newCoK(ch,:,2)))))%slice nr.2
    title(['ch no ' num2str(ch)]);
%     hold on
%     plot(tPulse, detrendedData/max(detrendedData),'g');
%     hold on,
%     plot(tResp,DetrendedDataTotal{2}/max(DetrendedDataTotal{2}),'r');
 end

figure(5)
for ch=24:31
    subplot(3,2,ch-23)
    plot(t,detrend(squeeze(newCoK(ch,:,2)))/max(detrend(squeeze(newCoK(ch,:,2)))))%slice nr.2
    title(['ch no ' num2str(ch)]);
%     hold on
%     plot(tPulse, detrendedData/max(detrendedData),'g');
%     hold on,
%     plot(tResp,DetrendedDataTotal{2}/max(DetrendedDataTotal{2}),'r');
end

L=length(newCoK);
nfft=2^nextpow2(L);
Fs=1/0.3;

figure(1)
for chan=1:6
    subplot(3,2,chan)
    ifftCoK = ifft(fft(squeeze(detrend(newCoK(chan,:,2))), nfft)/L);
    f = Fs/2*linspace(0,1,nfft/2+1);
    plot(f, 2*abs(ifftCoK(1:nfft/2+1)));
    title(['ch no ' num2str(chan)]);
end  

figure(2)
for chan=7:12
    subplot(3,2,chan-6)
    ifftCoK = ifft(fft(squeeze(detrend(newCoK(chan,:,2))), nfft)/L);
    f = Fs/2*linspace(0,1,nfft/2+1);
    plot(f, 2*abs(ifftCoK(1:nfft/2+1)));
    title(['ch no ' num2str(chan)]);
end  

L=length(newCoK);
nfft=2^nextpow2(L);
Fs=1/0.3;

%% Fourier analysis of each chanel using the periodogram method

figure (1)
for chan=1:6
    subplot(3,2,chan)
    Hs=spectrum.periodogram;      % Use default values
    psd(Hs,squeeze(detrend(newCoK(chan,:,3))),'Fs',Fs)
    title(['Periodogram ch no ' num2str(chan)]);
end  

figure(2)
for chan=7:12
    subplot(3,2,chan-6)
    Hs=spectrum.periodogram;      % Use default values
    psd(Hs,squeeze(detrend(newCoK(chan,:,3))),'Fs',Fs)
    title(['Periodogram ch no ' num2str(chan)]);
end  

figure(3)
for chan=13:19
    subplot(3,2,chan-12)
    Hs=spectrum.periodogram;      % Use default values
    psd(Hs,squeeze(detrend(newCoK(chan,:,3))),'Fs',Fs)
    title(['Periodogram ch no ' num2str(chan)]);
end  

figure(4)
for chan=20:25
    subplot(3,2,chan-19)
    Hs=spectrum.periodogram;      % Use default values
    psd(Hs,squeeze(detrend(newCoK(chan,:,3))),'Fs',Fs)
    title(['Periodogram ch no ' num2str(chan)]);
end  

figure(5)
for chan=26:31
    subplot(3,2,chan-25)
    Hs=spectrum.periodogram;      % Use default values
    psd(Hs,squeeze(detrend(newCoK(chan,:,3))),'Fs',Fs)
    title(['Periodogram ch no ' num2str(chan)]);
end 

%% 
figure
    subplot(3,2,1)
    ifftCoK = ifft(fft(squeeze(detrend(newCoK(1,:,2))), nfft)/L);
    f = Fs/2*linspace(0,1,nfft/2+1);
    plot(f, 2*abs(ifftCoK(1:nfft/2+1)));
    title(['ch no ' num2str(1)]);
    
      subplot(3,2,2)
    ifftCoK = ifft(fft(squeeze(detrend(newCoK(2,:,2))), nfft)/L);
    f = Fs/2*linspace(0,1,nfft/2+1);
    plot(f, 2*abs(ifftCoK(1:nfft/2+1)));
    title(['ch no ' num2str(2)]);
    
      subplot(3,2,3)
    ifftCoK = ifft(fft(squeeze(detrend(newCoK(6,:,2))), nfft)/L);
    f = Fs/2*linspace(0,1,nfft/2+1);
    plot(f, 2*abs(ifftCoK(1:nfft/2+1)));
    title(['ch no ' num2str(6)]);
    
      subplot(3,2,4)
    ifftCoK = ifft(fft(squeeze(detrend(newCoK(7,:,2))), nfft)/L);
    f = Fs/2*linspace(0,1,nfft/2+1);
    plot(f, 2*abs(ifftCoK(1:nfft/2+1)));
    title(['ch no ' num2str(7)]);
    
      subplot(3,2,5)
    ifftCoK = ifft(fft(squeeze(detrend(newCoK(17,:,2))), nfft)/L);
    f = Fs/2*linspace(0,1,nfft/2+1);
    plot(f, 2*abs(ifftCoK(1:nfft/2+1)));
    title(['ch no ' num2str(17)]);
    
      subplot(3,2,6)
    ifftCoK = ifft(fft(squeeze(detrend(newCoK(5,:,2))), nfft)/L);
    f = Fs/2*linspace(0,1,nfft/2+1);
    plot(f, 2*abs(ifftCoK(1:nfft/2+1)));
    title(['ch no ' num2str(5)]);

%% Fourier analysis each channel individually
figure(1)
for chan=1:6
    subplot(3,2,chan)
    ifftCoK = ifft(fft(squeeze(detrend(newCoK(chan,:,2))), nfft)/L);
    f = Fs/2*linspace(0,1,nfft/2+1);
    plot(f, 2*abs(ifftCoK(1:nfft/2+1)));
    title(['ch no ' num2str(chan)]);
end  

figure(2)
for chan=7:12
    subplot(3,2,chan-6)
    ifftCoK = ifft(fft(squeeze(detrend(newCoK(chan,:,2))), nfft)/L);
    f = Fs/2*linspace(0,1,nfft/2+1);
    plot(f, 2*abs(ifftCoK(1:nfft/2+1)));
    title(['ch no ' num2str(chan)]);
end  

figure(3)
for chan=13:19
    subplot(3,2,chan-12)
    ifftCoK = ifft(fft(squeeze(detrend(newCoK(chan,:,2))), nfft)/L);
    f = Fs/2*linspace(0,1,nfft/2+1);
    plot(f, 2*abs(ifftCoK(1:nfft/2+1)));
    title(['ch no ' num2str(chan)]);
end  

figure(4)
for chan=20:25
    subplot(3,2,chan-19)
    ifftCoK = ifft(fft(squeeze(detrend(newCoK(chan,:,2))), nfft)/L);
    f = Fs/2*linspace(0,1,nfft/2+1);
    plot(f, 2*abs(ifftCoK(1:nfft/2+1)));
    title(['ch no ' num2str(chan)]);
end  

figure(5)
for chan=26:31
    subplot(3,2,chan-25)
    ifftCoK = ifft(fft(squeeze(detrend(newCoK(chan,:,2))), nfft)/L);
    f = Fs/2*linspace(0,1,nfft/2+1);
    plot(f, 2*abs(ifftCoK(1:nfft/2+1)));
    title(['ch no ' num2str(chan)]);
end 



