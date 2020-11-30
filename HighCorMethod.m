%% HighCor method: 


function [rHighCor, CorHighCor, newMAGHighCor, CorrMAP, SelectedVoxels]=HighCorMethod(MAG,dMAG,PH,uPH,tSTDmask_MAG)


%%3. Pearson Correlation between the phase and the magnitude images on the voxels with the higher tSTD:

%select the voxels with higher tSTD:
[row, col]=find(tSTDmask_MAG);
corr_index=[row col];

CorrMAP=zeros(size(PH,1),size(PH,2));
SelectedVoxels=zeros(size(PH,1),size(PH,2));
t=0;
for i=1:size(PH,1)
    for j=1:size(PH,2)
        %select only the voxels that were selected with the CompCor method
        if find(ismember(corr_index,[i j],'rows')==true)
            [corr]=corrcoef(squeeze(uPH(i,j,7,:)),squeeze(dMAG(i,j,7,:)));
            CorrMAP(i,j)=corr(1,2);
            if abs(CorrMAP(i,j))>0.4;
                SelectedVoxels(i,j)=CorrMAP(i,j);
                t=t+1;
                timeSeries(t,:)=squeeze(dMAG(i,j,7,:));
            else 
                SelectedVoxels(i,j)=0;
            end 
        else 
            CorrMAP(i,j)=0;
        end  
    end 
end

% Applying SVD on the timeSeries of the chosen voxels

[~,SVDMatrix,V]=svd(timeSeries(:,:));
lambda=diag(SVDMatrix);
variances=lambda.*lambda;
PCs=6;
VV=V(:,1:PCs);
Y=VV'*timeSeries';
XX=VV*Y;
XX=XX';
rHighCor=squeeze(XX(1,:))';

%% Regress out the PCA of the timeSeries of the chosen voxels

CorHighCor=zeros(size(MAG,1),size(MAG,2),size(MAG,4));
newMAGHighCor=zeros(size(MAG,1),size(MAG,2), size(MAG,4));

for i=1:size(MAG,1)
    for j=1:size(MAG,2)
        CorHighCor(i,j,:)=regress(squeeze(dMAG(i,j,7,:)),rHighCor);
        newMAGHighCor(i,j,:)=mean(squeeze(dMAG(i,j,7,:)))+(squeeze(dMAG(i,j,7,:))-squeeze(CorHighCor(i,j,:)));
    end 
end

%% figures to check

figure
bar(variances(1:30));

figure
for i=1:6
    plot(squeeze(XX(i,:))');
    hold on
end 

figure
subplot(131)
imagesc(rot90(tSTDmask_MAG(:,:)))
axis off
subplot(132)
imagesc(rot90(CorrMAP(:,:)))
axis off
subplot(133)
imagesc(rot90(SelectedVoxels(:,:)))
axis off
end
