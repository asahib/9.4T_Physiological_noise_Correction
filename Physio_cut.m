%% cut the physio files 

function [new_Physio]=sync_Physio(Physio, dcm_timings)

%1. Get the Acquisition start time of the first time point and the last time point of the volume (dcm header)

hdr1=spm_dicom_headers(spm_select);
hdr2=spm_dicom_headers(spm_select);
TimePoint1_AcqTime=hdr1{1,1}.AcquisitionTime;
TimePointlast_AcqTime=hdr2{1,1}.AcquisitionTime;

%2. Cut the beginning of the resp and pulse checking the LogStartMDHTime for resp and pulse

RespStart=Physio

end
