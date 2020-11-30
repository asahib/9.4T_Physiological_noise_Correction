spmdir='/p/joana';
subjects={'MSBOLD03' 'MSBOLD04' 'MSBOLD05' 'MSBOLD12' };
%subjects={'MSBOLD03'};
subdirs={'session1';'session2'};
%subdirs={'session1' };

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
                hdr2=spm_dicom_headers(spm_select('FPList', dcmdir,'.*\.0301.*\.IMA'));
                
                TimePoint1_AcqTime=hdr1{1,1}.AcquisitionTime;
                TimePointlast_AcqTime=hdr2{1,1}.AcquisitionTime;
                
                dcm_timings(subj,subdir,1)=TimePoint1_AcqTime;
                dcm_timings(subj,subdir,2)=TimePointlast_AcqTime;
                
                dcm_timings_subj=dcm_timings(subj,subdir,:);
                
                %select physiological data files
                Resp_raw=read_physio(spm_select('FPList',physiodir, '^resptxt_.*\.resp'));
                Pulse_raw=read_physio(spm_select('FPList', physiodir,'^pulsetxt_.*\.puls'));
                Physio_raw={Resp_raw,Pulse_raw};
                
                %cut the physio files to match the dcms
                [Physio]=sync_Physio(Physio_raw,dcm_timings(subj,subdir,:),TR,dtr);
                
                %Interpolation of respiration signal
                %[respiration_model,onset]=RecreateRespNew(squeeze(Physio(1,:)),dtr,[]);
                
                %determine RETROICOR regressor
                rRetroicor=RetroicorReg(Physio,600);
                
                %put NaNs to zeros
                rRetroicor(find(isnan(rRetroicor)))=0;
                
                %detrend RETROICOR regressor
                for i=1:size(rRetroicor,2)
                    R(:,i)=detrend(rRetroicor(:,i),'constant');
                end
                
                %determine covariates matrix
                 save([Covdir '/' 'rRetroicor_' subjects{subj} '_' subdirs{subdir} '.mat'],'R')
                 save([Covdir '/' 'dcmtimings_' subjects{subj} '_' subdirs{subdir} '.mat'],'dcm_timings_subj')
                 
                 clear rRetroicor R Resp_raw Pulse_raw Physio_raw Physio hdr1 hdr2
                 
    end
end
