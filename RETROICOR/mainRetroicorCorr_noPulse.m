spmdir='/p/joana';
subjects={'MSBOLD02' 'MSBOLD10' 'MSBOLD11'};
%subdirs={'session1';'session2'};
subdirs={'session1' 'session2'};

TR=1;
dtr=49.8;
for subj=1:length(subjects)
    for subdir=1:length(subdirs)
                        
                %data directory
                dcmdir=fullfile(spmdir,'DCM_MSBOLD',subjects{subj}, subdirs{subdir});
                Covdir=fullfile(spmdir,'NIFTI_MSBOLD',subjects{subj},'Covariates', subdirs{subdir});
                physiodir=fullfile(spmdir,'NIFTI_MSBOLD',subjects{subj},'physio',subdirs{subdir});
                
                %select the first and last dcm time points for each session per subject
                hdr1=spm_dicom_headers(spm_select('FPList', dcmdir,'.*\.0002.*\.IMA'));
                hdr2=spm_dicom_headers(spm_select('FPList', dcmdir,'.*\.0601.*\.IMA'));
                
                TimePoint1_AcqTime=hdr1{1,1}.AcquisitionTime;
                TimePointlast_AcqTime=hdr2{1,1}.AcquisitionTime;
                
                dcm_timings(subj,subdir,1)=TimePoint1_AcqTime;
                dcm_timings(subj,subdir,2)=TimePointlast_AcqTime;
                
                %select physiological data files
                Resp_raw=read_physio(spm_select('FPList',physiodir, '^resptxt_.*\.resp'));
                
                %cut the physio files to match the dcms
                [Physio]=sync_Physio_noPulse(Resp_raw,dcm_timings(subj,subdir,:),TR,dtr);
                
                %Interpolation of respiration signal
                %[respiration_model,onset]=RecreateRespNew(squeeze(Physio(1,:)),dtr,[]);
                
                %determine RETROICOR regressor
                rRetroicor=RetroicorReg_noPulse(Physio);
                
                %put NaNs to zeros
                rRetroicor(find(isnan(rRetroicor)))=0;
                
                %detrend RETROICOR regressor
                for i=1:size(rRetroicor,2)
                    R(:,i)=detrend(rRetroicor(:,i),'constant');
                end
                
                %determine covariates matrix
                 save([Covdir '/' 'rRetroicor_' subjects{subj} '_' subdirs{subdir} '.mat'],'R')
                 
                 clear rRetroicor R Resp_raw Pulse_raw Physio_raw Physio hdr1 hdr2
                 
    end
end
