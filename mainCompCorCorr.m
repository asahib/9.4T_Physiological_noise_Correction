spmdir='/p/joana/MSBOLD_paper';
%subjects={'sub1'};
subjects={'sub1' 'sub2' 'sub3' 'sub4' 'sub5' 'sub6' 'sub7' 'sub8' 'sub9'};
%subjects={'MSBOLD02'};
subdirs={'session1' 'session2'};

TR=1;
dtr=49.8;
for subj=1:length(subjects)
    for subdir=1:length(subdirs)
        
        %data directory
        voldir=fullfile(spmdir,subjects{subj},'Functionals','Preprocessed_images', subdirs{subdir},'PSF');
        Covdir=fullfile(spmdir,'Regressors', 'CompCor', '0_06PC_6comp', subdirs{subdir},'PSF');
        
        %select EPI Volume
        Vol=spm_read_vols(spm_vol(spm_select('FPList', voldir,'^rf.*\.nii')));
        
        %Compute tSTD Mask (voxels with 2% highest tSTD)
        mask=tCompCor_mask(Vol,0.00009);
        
        %get the timeseries of the voxels from the mask
        mask_tseries=get_tseries(Vol,mask);
        mask_tseries(all(mask_tseries==0,2),:)=[];
        
        %detrend the timeseries
        for i=1:size(mask_tseries,1)
            mask_tseries_det(i,:)=detrend(mask_tseries(i,:),'constant');
        end
        
        %Compute PCA analysis using SVD
        [u s v] = svd(mask_tseries_det');
        % get the first 6 eigenvectors to use as regressors
        R=u(:,1:5);
        
        %sve covariates matrix
        save([Covdir '/' 'rCompCor_0_001PC_5comp_'  subjects{subj} '_' subdirs{subdir} '_PSF.mat'],'R')
        
        %save tSTD mask used
        mask=struct('vol',mask);
        MRIwrite(mask,[Covdir '/' 'tSTD_mask_0_001PC' subjects{subj} '_' subdirs{subdir} '_PSF.nii' ] ,'double');
        
        clear R mask_tseries mask_tseries_det Vol mask u s v
    end
end
