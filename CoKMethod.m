%% CoK Method (information from the raw data of the images)

function [rCoK, CorCoK, newMAGCoK]=CoKMethod(MAG,dMAG,rawData)

[CoK, newCoK, image_obj] =GetCoK_modified_Gisela('meas_MID1078_mzBOLD_1190_4sl_0_9mSineNoDiCoTE20P4TR300_FID1703.dat');
dummy = squeeze(image_obj{210,:,54,1,:});

%% Clean the signal correspondent to the center of k space 

%% Reduce data dimensions using SVD


%% Regress the signal of one chanel out

CoKStab050012_ch20=zeros(size(VStab050012,1),size(VStab050012,2));
 newVStab050012CoKch20Corrected=zeros(size(VStab050012,1),size(VStab050012,2), 1000);
for i=1:size(VStab050012,1)
    for j=1:size(VStab050012,2)
        CoKStab050012_ch20(i,j)=regress(squeeze(VStab050012(i,j,3,1:1000)),squeeze(newCoK(20,:,3))');
        newVStab050012CoKch20Corrected(i,j,:)=mean(squeeze(VStab050012(i,j,3,1:1000)))+(squeeze(VStab050012(i,j,3,1:1000))-(CoKStab050012_ch20(i,j)*squeeze(newCoK(20,:,3)))');

    end 
end

end 