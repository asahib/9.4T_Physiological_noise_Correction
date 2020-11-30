function [respiration_model,onset]=RecreateRespNew(dat,dtr,mxResp)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%function [respiration_model,onset]=RecreateRespNew(dat,dtr,mxResp);
%recreate respiration data

% function [respiration_model,onset]=RecreateRespNew(dat,mxResp,dtr)
% dat: respiration trace
% dtr: sampling frequency pr second (typically 50)
% mxResp: maxima in respiration trace found by scanner

% this function first fits a spline baseline, which is subtracted from the
% input. The new data is then rescaled to the range 0-100, so that the
% amplitude threshold for the peak fitting works OK
%
% respiration data is the recreated based on iterative peak fitting and
% GLM fitting to find the correct depth of breath also in presence of
% clipping
%
%Gisela Hagberg, High Field Magneti Resonance, MPI for biological
%Cybernetcis and University Hospital Tuebingen
%Marta Bianciardi, Martino Center, MGH Boston
%
%while at Santa Lucia Foundation Rome, 2005-2010
%calls: 
%FindFirstTemplate
%GenRegrModelTuk
%calcola_correlation_ECG
Fc=1/dtr;
%find peaks
time=[0:length(dat)-1]*dtr;
if isempty(mxResp)
  [pks mxResp]=findpeaks(dat,'minpeakdistance',45,'minpeakheight',2048);

end
plot(time,dat,'-',time(mxResp),dat(mxResp),'or')  
fprintf('are all maxima found then press retur, else interrupt and modify mxResp')
pause


%find knots for baseline correction, Method I : timepoints for inspiration
%discovered by scanner

%for k=1:length(mxResp)-1
%    [Y,I]=min(dat(mxResp(k):mxResp(k+1)));
%    idx(k)=mxResp(k)+I;
%end;
%dknots=idx(dat(idx)>(mean(dat(idx))-1.96*std(dat(idx))));

%find knots for baseline correction, Method II : automatic

[pks dknots]=findpeaks(-1*dat,'minpeakdistance',25,'minpeakheight',-0.25*max(dat));
idx=find(dat==0);
idx=idx(find(diff([0 idx])>1));
dknots=sort([dknots idx]);
dknots=dknots(find(diff([0 dknots])>1));


y=dat(dknots);
my=mean(y)
dknots=[1 dknots length(dat)];
y=[my y my];

 baseline=interp1(dknots,y,[1:length(dat)]);
newdat=dat-baseline;
factor=100/max(newdat);
newdat=newdat*factor;
baseline=baseline*factor;

thresh.grad=50;%no of timepoints for peakfinding mean(diff(mxResp))*dtr;
thresh.ampli=20;%20% of max

fprintf('1. zoom in to one respiration cycle; 2. press return; 3. select the two adjacent minima of one cycle');



[newRESP_template, onset_newRESP]=FindFirstTemplateResp(newdat,thresh,Fc)

%generates model, no detrend of template, corrects possible local problems
%in baseline correction
[regressor_model_conv, regressor_model, beta]=GenRegrModelTuk(newdat',onset_newRESP,newRESP_template,[]);

%now let's replace missing respiration information
onset=onset_newRESP;
onsbeta=beta;

%finally generate the new respiration trace!
tukey = tukeywin(size(newRESP_template,2), 0.05)';
stick = zeros(1,length(newdat));
stick(onset) = onsbeta;
respiration_model = conv(stick, newRESP_template.*tukey);    
respiration_model=respiration_model(1:length(dat))/factor+baseline;
return;




