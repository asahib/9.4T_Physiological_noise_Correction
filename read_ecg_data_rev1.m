function [ecg, trigger, stats]=read_ecg_data_rev1(filename, plotflag);
% function [ecg, trigger, stats]=read_ecg_data_rev1(filename, plotflag);
%
% function to read ecg data from Siemens PMU

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

% use dialog
idir='.'; % starting directory for UI dialog
if nargin==0
    [file,path]=uigetfile([idir,'*.ecg'],'Select ECG file:');
    if isequal(file,0) | isequal(path,0)
       return;
    end
    filename=[path,file];
end

% default plot flag to be off
if nargin<2;
    plotflag=1;
end

[p, ecg_filename]= fileparts(filename);

fid = fopen(filename);
i=0;
while ~feof(fid) 
   line = fgetl(fid);
   i=i+1;
   if isempty(line), break, end
   data{i}=line;
end
fclose(fid);
d1=data{1};
ecg=d1;
index = findstr(ecg,'maxLimitAVF');
tmp=find(isspace(ecg(index:index+100)));
index=index+tmp(2);


ecg=str2num(ecg(index:end));
if nargout==3
    stats=strvcat(data{2:end});
end

% discard 1st 5 header values
% header=ecg(1:5);


% 2 1 2 40 280 628 12175 628 12175 628 12175 628 12175 628 12175 628 12175
% 628 5002 eTriggerMethod: 1, minLimitCh1: 39, maxLimitCh1: 398, minLimitAVF: 22, maxLimitAVF: 224

trigger=zeros(size(ecg));
start_triggers=find(ecg==5000);
for i=1:length(start_triggers)
    if start_triggers(i)==1
        trigger(start_triggers(i)+1)=1;
    else
        trigger(start_triggers(i)-1)=1;
    end
end

trigger=trigger(~or(ecg==5000,ecg==6000));

ecg=ecg(~or(ecg==5000,ecg==6000));


% check if 1st value is channel 1 or 2
if ecg(1) > 6000 % then must be channel 2
    ecg=ecg(2:end);
    trigger=trigger(2:end);
end
% check for even number of samples
if mod(length(ecg),2)~=0
    ecg=ecg(1:end-1);
    trigger=trigger(1:end-1);
end
% separate channels
ecg=reshape(ecg,[2 length(ecg)/2]);
trigger=reshape(trigger,[2 length(trigger)/2]);
trigger=or(trigger(1,:),trigger(2,:));

% remove bias
ecg(1,:)=ecg(1,:)-2048;
ecg(2,:)=ecg(2,:)-10240;

% delete last sample which is sometimes bad (??)
ecg=ecg(:,1:end-1);
trigger=trigger(:,1:end-1);

if plotflag==1
%     figure
%     plot(ecg(1,:),ecg(2,:));
%     hold on
%     plot(ecg(1,trigger),ecg(2,trigger),'ro')
%     axis equal
    
    t=[1:size(ecg,2)]*0.0025; % 400 samples/second
    
    ymax=2048;
    h_fig=figure;
    h_axis(1)=subplot(2,1,1); plot(t,ecg(1,:));shg
    ecg_filename=strrep(ecg_filename,'_',' ');
    title(ecg_filename);
    hold on;plot(t(trigger),ecg(1,trigger),'ro');shg
%     axis([t(1) t(end) -ymax ymax]);
    h_axis(2)=subplot(2,1,2); plot(t,ecg(2,:));shg
    hold on;plot(t(trigger),ecg(2,trigger),'ro');shg
%     axis([t(1) t(end) -ymax ymax]);
    linkaxes(h_axis,'x');
    zoom xon
    pan(h_fig,'xon');
end

return

% 
% 1. lPmuECGModePub                       // ECG control mode 
%         METHOD_NONE        = 0x01 
%         METHOD_TRIGGERING  = 0x02 
%         METHOD_GATING      = 0x04 
%         METHOD_RETROGATING = 0x08 
%         METHOD_SOPE        = 0x10 
%         METHOD_ALL         = 0x1E 
% 2. lPmuADPub                            // arrhythmia detection mode 
%         don't care 
% 3. iPmuHighPrioTriggerSignal            // high prio trigger signal 
%         SIGNAL_NONE        = 0x01, 
%         SIGNAL_EKG         = 0x02, 
%         SIGNAL_PULSE       = 0x04, 
%         SIGNAL_EXT         = 0x08, 
%         SIGNAL_RESPIRATION = 0x10, 
%         SIGNAL_ALL         = 0x1E, 
%         SIGNAL_EKG_AVF     = 0x20 
% 4. ulECGGateOnCountPub          // time stamp counter gate ON ECG 
%                                         // ON:  100 ms  / 2.5 = 40  
% 5. ulECGGateOffCountPub         // time stamp counter gate OFF ECG 
%                                         // OFF: 700 ms  / 2.5 = 280 
%                                  
% NOTE:
% If the latest Siemens VCG WIP (Improvements S073) has been installed on
% scanner, the outputs stats will include start and stop time for MDH and
% MPCU to allow better synchronization with the raw files
%
% LogStartMDHTime:                
% LogStopMDHTime:                 
% LogStartMPCUTime:               
% LogStopMPCUTime:   

% 
% below is a description of the header and footer data for the VB13A:
%  
% description of the addtional header data between 5002 and 6002:
%  
% eTriggerMethod: 1, minLimitCh1: 17, maxLimitCh1: 176, minLimitAVF: 3, maxLimitAVF: 38
%  
% eTriggerMethod: 1 means "standard VCG" trigger method, 2 means "channel I only", 3 means
%  "channel AVF only".
%  The remaining 4 values are just for internal debugging purposes, and can be ignored.
% description of the footer between 5003 and 6003:
% 
% short term statistics (current short term averages at the end of the data logging):
% 
% ECG Freq Per: 0 0
% PULS Freq Per: 0 0
% RESP Freq Per: 0 0
% EXT Freq Per: 0 0
% 
% long term statistics (long term averages etc. from the begin of the last protocol till the end of the data logging):
% 
% ECG Min Max Avg StdDiff: 0 0 0 0
% PULS Min Max Avg StdDiff: 610 2503 959 19
% RESP Min Max Avg StdDiff: 0 0 0 0
% EXT Min Max Avg StdDiff: 0 0 0 0
% NrTrig NrMP NrArr AcqWin: 49 0 0 921
% 
% the following values allow for the alignment between the MR raw data and the logged signal values - times are given in msec
% 
% LogStartMDHTime: 56840670 // MDH timestamp at the start time of the signal logging
% LogStopMDHTime: 57049812 // MDH timestamp at the stop time of the signal logging
% LogStartMPCUTime: 56828815 // MPCU timestamp at the start time of the signal logging 
% LogStopMPCUTime: 57049935 // MPCU timestamp at the stop time of the signal logging

                                        