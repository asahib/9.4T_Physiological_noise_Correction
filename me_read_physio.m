function [res, meta]=me_read_physio(action, physiofile, flag)
% function [res, meta]=me_read_physio(action, physiofile, flag)
% read physio logfiles from Siemens
% output:
%           res                   time course of signal(s) nx1 or nx2(ecg)
%           meta                  meta informations
%           meta.filename         name of physiofile
%           meta.sampling         sampling intervall in us (microseconds)
%           meta.header           header infos <unknown>
%           meta.strings          cell array of cellstrings with info srings (5002 - 6002)
%           meta.stat             cell array of cellstrings with statistic srings (5003 - 6003)
%           meta.time             sample timepoints (in secs)
%           meta.time_e1          time of events 1 (5000)
%           meta.time_e2          time of events 2 (6000)
%           meta.time_strings     time of strings (5002 - 6002)
%           meta.trigger          time of Trigger events ("Trig:")
%           meta.trigger_off      relative offset time of last detected event (5000?)
%
%           meta.LogStartMDHTime  start of log time MDH  (msecs since midnight)
%           meta.LogStopMDHTime   stop of log time MDH   (msecs since midnight)
%           meta.LogStartMPCUTime start of log time MPCU (msecs since midnight)
%           meta.LogStopMPCUTime  stop of log time MPCU  (msecs since midnight)

if ~exist('action','var')
    action='old';
    % action='datavar';
end
if exist(action,'file')
    if exist('physiofile','var'), flag=physiofile; end
    physiofile=action;
    action='data';
end
if ~exist('physiofile','var')
    fprintf('No files specifiesd!\n');
    return
end
if ~exist('flag','var')
    flag=0;
    %flag=1;  % plot
    %flag=2; % print out
end
meta.filename=physiofile;
switch lower(action)
    case 'data'
        [pt, nm, ext]=fileparts(physiofile);
        switch lower(ext)
            case '.ecg'
                samplingrate=400; % Hz
                nchan=4;
            case '.ext'
                samplingrate=200; % Hz
            otherwise
                samplingrate= 50; % Hz
        end
        meta.sampling = 1000000/samplingrate;    % sampling intervall in us (microseconds)
        %fid=fopen(fullfile(physiopath,subdirs{1},dir_info(1).name),'r');
        fid=fopen(physiofile,'r');
        data=textscan(fid,'%s');
        fclose(fid);
        data=data{1};
        mask=ones(length(data),1);
        tmask=zeros(size(mask));
        % events 1
        event_index1  = find(strcmp(data,'5000'));
        mask(event_index1)= 0;
        % events 2
        event_index2  = find(strcmp(data,'6000'));
        mask(event_index2)= 0;
        % events 3
        event_index3  = find(strcmp(data,'7000'));   % not shure if any exist
        mask(event_index3)= 0;
        start_str_index  = find(strcmp(data,'5002'));
        stop_str_index   = find(strcmp(data,'6002'));
        start_stat_index = find(strcmp(data,'5003'));
        if isempty(start_stat_index), start_stat_index = find(strcmp(data,['5003' char(13)])); end
        stop_stat_index  = find(strcmp(data,['6003']));
        % header info <unknown>
        nheader=4;
        % meta.header=str2num(char(data(1:start_str_index(1)-1)));
        %mask(1:start_str_index(1))=0;
        meta.header=str2num(char(data(1:nheader)));
        mask(1:nheader)=0;
        % info strings
        trig_index=[];
        trig_offset=[];
        for s=1:length(start_str_index)
            meta.strings{s}=data(start_str_index(s)+1:stop_str_index(s)-1);
            mask(start_str_index(s):stop_str_index(s))=0;
            if strcmp(meta.strings{s}{1},'Trig:')
                trig_index = [trig_index; start_str_index(s)];
                trig_offset = [trig_offset; str2double(meta.strings{s}{3})];
            end; 
            if strcmp(meta.strings{s}{1},'1SecToStart')
                meta.OneSecToStart=s; 
            end;
            if strcmp(meta.strings{s}{1},'2SecToStart')
                meta.TwoSecToStart=s; 
            end;
            if strcmp(meta.strings{s}{1},'3SecToStart')
                meta.ThreeSecToStart=s; 
            end;
        end;
        % determine number of channels (ecg)
        if strcmp(meta.strings{1}{1},'LOGVERSION')
            switch str2double(meta.strings{1}{2})
                case 2
                    nchan=2;
                case 3
                    nchan=4;
            end
        end
        % statistic strings
        for s=1:length(start_stat_index)
            meta.stat{s}=data(start_stat_index(s)+1:stop_stat_index(s)-1);
            mask(start_stat_index(s):stop_stat_index(s))=0;
        end;
        % residual data values
        res=str2num(char(data(mask==1)));
        % time points
        tmask(mask==1)=1:length(find(mask==1));
        %if exist('vect','var')
        if ~isempty(res)
            if samplingrate==400    % ecg file
                switch nchan
                    case 2
                        if  mod(length(res),2)      % odd number of samples
                            res=[res; res(end-1)];
                            tmpind=find(mask>0);
                            tmask=[tmask(1:tmpind(end)); tmask(tmpind(end):end)];
                        end
                        tmask=fix((tmask+1)/2);     % correct sample number
                        res=reshape(res,2,length(res)/2)';
                        res(:,1)=res(:,1)-2048;
                        res(:,2)=res(:,2)-10240;
                    case 3
                        if  mod(length(res),3)      % odd number of samples
                            res=[res; res(end-1)];
                            tmpind=find(mask>0);
                            tmask=[tmask(1:tmpind(end)); tmask(tmpind(end):end)];
                        end
                        tmask=fix((tmask+1)/3);     % correct sample number
                        res=reshape(res,3,length(res)/3)';
                        res(:,1)=res(:,1)-2048;
                        res(:,2)=res(:,2)-10240;
                        res(:,3)=res(:,3)-18432;
                    case 4
                        if  mod(length(res),4)      % odd number of samples
                            res=[res; res(end-1)];
                            tmpind=find(mask>0);
                            tmask=[tmask(1:tmpind(end)); tmask(tmpind(end):end)];
                        end
                        tmask=fix((tmask+1)/4);     % correct sample number
                        res=reshape(res,4,length(res)/4)';
                        res(:,1)=res(:,1)-2048;  %  2 * 1024  [    0 -  4096]
                        res(:,2)=res(:,2)-10240; % 10 * 1024  [ 8192 - 12288]
                        res(:,3)=res(:,3)-18432; % 18 * 1024  [16384 - 20480]
                        res(:,4)=res(:,4)-26624; % 26 * 1024  [24576 - 28672]
                end
            end
            % time 
            meta.time=((1:length(res))/samplingrate)';
            %meta.time=((0:length(res)-1)/samplingrate)';
            if ~isempty(event_index1), meta.time_e1=tmask(event_index1-1)/samplingrate; end
            if ~isempty(event_index2), meta.time_e2=tmask(event_index2-1)/samplingrate; end
            meta.time_strings=tmask(start_str_index-1)/samplingrate;
            if isfield(meta,'OneSecToStart'),   meta.time_OneSecToStart   = meta.time_strings(meta.OneSecToStart);   end
            if isfield(meta,'TwoSecToStart'),   meta.time_OneSecToStart   = meta.time_strings(meta.TwoSecToStart);   end
            if isfield(meta,'ThreeSecToStart'), meta.time_ThreeSecToStart = meta.time_strings(meta.ThreeSecToStart); end
            % meta.time_statistic=1000*tmask(start_stat_index-1)/samplingrate;
            if ~isempty(trig_index), meta.trigger=tmask(trig_index-1)/samplingrate; end
            if ~isempty(trig_offset), meta.trigger_off=trig_offset/samplingrate; end
            % log time (msecs since midnight)-> secs
            meta.LogStartMDHTime  =  str2double(meta.stat{1}{end-6})/1000;
            meta.LogStopMDHTime   =  str2double(meta.stat{1}{end-4})/1000;
            meta.LogStartMPCUTime =  str2double(meta.stat{1}{end-2})/1000;
            meta.LogStopMPCUTime  =  str2double(meta.stat{1}{end})/1000;
            meta.time2=((1:length(res))*(meta.LogStopMDHTime-meta.LogStartMDHTime)/(length(res)-1))';
            %meta.time2=((0:length(res)-1)*(meta.LogStopMDHTime-meta.LogStartMDHTime)/length(res))';
            %meta.sampling = (meta.LogStopMDHTime-meta.LogStartMDHTime)/length(res);    % sampling intervall in us (microseconds)
            if flag
                figure
                %plot(meta.time, res)
                plot(meta.time2, res)
                hold on
                if isfield(meta,'trigger'), plot(meta.trigger, ones(size(meta.trigger))*1024,'r+'), end
                if isfield(meta,'time_e1'), plot(meta.time_e1, ones(size(meta.time_e1))*850,'b+'), end
                if isfield(meta,'time_e2'), plot(meta.time_e2, ones(size(meta.time_e2))*850,'g+'), end
                if isfield(meta,'OneSecToStart'),   plot(meta.time_OneSecToStart*[1 1],     [0 4095],'r'), end
                if isfield(meta,'OneSecToStart'),   plot(meta.time_OneSecToStart*[1 1]+1,   [0 4095],'m'), end
                if isfield(meta,'TwoSecToStart'),   plot(meta.time_TwoSecToStart*[1 1],     [0 4095],'r'), end
                if isfield(meta,'TwoSecToStart'),   plot(meta.time_TwoSecToStart*[1 1]+1,   [0 4095],'m'), end
                if isfield(meta,'ThreeSecToStart'), plot(meta.time_ThreeSecToStart*[1 1],   [0 4095],'r'), end
                if isfield(meta,'ThreeSecToStart'), plot(meta.time_ThreeSecToStart*[1 1]+3, [0 4095],'m'), end
                hold off
                title(physiofile,'interpreter','none');
                drawnow;
            end
            if flag > 1
                fprintf('\nMDH    Time: %8.3f s\n', (meta.LogStopMDHTime-meta.LogStartMDHTime));
                fprintf('MPCU   Time: %8.3f s\n', (meta.LogStopMPCUTime-meta.LogStartMPCUTime));
                fprintf('Sample Time: %8.3f s\n', length(res)/samplingrate);
                fprintf('-> Sample Rate (MDH):  %6.3f Hz (%6.3f ms)\n', length(res)/(meta.LogStopMDHTime-meta.LogStartMDHTime),  1000*(meta.LogStopMDHTime-meta.LogStartMDHTime)/length(res));
                fprintf('-> Sample Rate (MPCU): %6.3f Hz (%6.3f ms)\n', length(res)/(meta.LogStopMPCUTime-meta.LogStartMPCUTime),  1000*(meta.LogStopMPCUTime-meta.LogStartMPCUTime)/length(res));
                if isfield(meta,'trigger'), fprintf('TR = %6.3f s\n', mean(diff(meta.trigger))/1000); end
                fprintf('Start Time: %s\n', datestr(meta.LogStartMDHTime / (24*60*60),'HH:MM:SS.FFF'));
                fprintf('Stop  Time: %s\n', datestr(meta.LogStopMDHTime  / (24*60*60),'HH:MM:SS.FFF'));
            end
        else
            res=-1;
        end
    case 'datavar'
        %res=var(me_physio_neu('data', physiofile));
        res=var(me_read_physio('data', physiofile));
    otherwise
        fprintf('Unknown action %s\n',action);
            res=-1;
end


% Tr =1700ms, 34 slice 3x3x4+1
