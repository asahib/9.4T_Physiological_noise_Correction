%% cut the physio files 

function [Physio ]=sync_Physio_noPulse(resp_raw,dcm_timings,TR,dtr)

% resp_raw=Physio_raw{1}; %put respiration file
% pulse_raw=Physio_raw{2}; %put cardiac file

RespStart=resp_raw.resp.LogStartMDHTime;
%PulseStart=pulse_raw.pulse.LogStartMDHTime;
RespStop=resp_raw.resp.LogStopMDHTime;
%PulseStop=pulse_raw.pulse.LogStopMDHTime;

dcmStart=dcm_timings(1); %in seconds
dcmLast=dcm_timings(2);

%1. Cut the beginning of the resp and pulse checking the LogStartMDHTime for resp and pulse

x=dcmStart-RespStart;%extra time in respiration before acq start
cutresp=ceil((x)*dtr);%nr. of points in the respiration vector have to be cut in the beginning

% x=dcmStart-PulseStart;%extra time in pulse before acq start
% cutpulse=ceil((x)*dtr);%nr. of points in the pulse vector have to be cut in the beginning

%2. total time of the dicom acquisition

dummy=dcmLast-dcmStart;
TotalTime_dcm=dummy;
nrPhysioPoints=floor((TotalTime_dcm)*dtr); %49.8 is the sampling frequency - it is tru for resp and pulse

%3. knowing 1. and 2. cut the vectors accordingly

resp=resp_raw.respiration_values(cutresp:cutresp+nrPhysioPoints);
%pulse=pulse_raw.pulse_values(cutpulse:cutpulse+nrPhysioPoints);

Physio(1,:)=resp';
%Physio(2,:)=pulse';


end
