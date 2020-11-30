%% NewCoK regression Method: Regressing the center of k-space from the phase and amplitude images (raw data not needed)

function [rNewCoK,rNewCoKdet, newMAGNewCoK]=NewCoKMethod(PH,uPH,MAG,dMAG)
%function [rNewCoK newMAGNewCoK]=NewCoKMethod(PH,uPH,MAG,dMAG)


%get the center of the k-space
rNewCoK=zeros(size(MAG,3),size(MAG,4));
ksp=zeros(size(MAG,3),size(MAG,1),size(MAG,2), size(MAG,4));
for time=1:size(MAG,4)
    for slice=1:size(MAG,3)
    MAGSlc=squeeze(MAG(:,:,slice,time));
    uPHSlc=squeeze(uPH(:,:,slice,time));
    CSlc=MAGSlc.*exp(-1i*uPHSlc);
    ksp(slice,:,:,time)=fftshift(ifft2(CSlc));
    rNewCoK(slice,time)=angle(squeeze(ksp(slice,106,107,time)));
    rNewCoK1(slice,time)=angle(squeeze(ksp(slice,106,106,time)));
    end  
end  

clear slice
%Detrend and clean rNewCoK to try to catch respiratiopnn and heart beat
rNewCoKdet=zeros(size(MAG,3),size(MAG,4));
for slice=1:size(MAG,3)
    rNewCoKdet(slice,:) = detrendnonlin(rNewCoK1(slice,:),3);
    rNewCoKdet(slice,:) = detrendnonlin(rNewCoKdet(slice,:),3);
end 
clear slice
%Regress out the center of k-space from the detrended magnitude image
b=zeros(size(MAG,3),2);
newMAGNewCoK=zeros(size(MAG,1),size(MAG,2),size(MAG,3),size(MAG,4));


% for slice=1:size(dMAG,3)
%     for i=1:size(dMAG,1)
%         for j=1:size(MAG,2)
%             b(slice,:)=pinv([rNewCoK1(slice,:); ones(1,size(MAG,4))])'*squeeze(dMAG(i,j,slice,:));
%             newMAGNewCoK(i,j,slice,:)=(squeeze(dMAG(i,j,slice,:))-b(slice,1)*rNewCoK1(slice,:)');
%         end  
%     end 
% end  
 dummy=zeros(size(MAG,1),size(MAG,2),size(MAG,3),size(MAG,4));
 
for slice=1:size(dMAG,3)
 X=[detrend(squeeze(rNewCoK1(slice,:))); ones(1,size(MAG,4))];
 pinvX=pinv(X);
    for i=1:size(dMAG,1)
        for j=1:size(MAG,2)
            b(slice,:)=pinvX'*squeeze(dMAG(i,j,slice,:));
            dummy(i,j,slice,:)=(squeeze(detrend(dMAG(i,j,slice,:)))-b(slice,1)*X(1,:));
            newMAGNewCoK(i,j,slice,:)=mean(dMAG(i,j,slice,:),4)+dummy(i,j,slice,:);
        end  
    end 
end  
clear slice

%plots to check

figure
subplot(121)
imagesc(rot90(mean(PH(:,:,7),4)));
axis off
title('Raw Phase images')
subplot(122)
imagesc(rot90(mean(uPH(:,:,7,:),4)));
axis off
title('Phase images unwraped in time')


figure(2)
for slice=1:size(MAG,3)
    subplot(4,4,slice)
    imagesc(abs(squeeze(ksp(slice,:,:,100))));
end 

clear slice
figure
for slice=1:size(MAG,3)
    plot(squeeze(rNewCoK1(slice,:)))
     hold on
    plot(squeeze(rNewCoKdet(slice,:)), 'r')
end 