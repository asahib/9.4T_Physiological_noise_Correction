function [regressor_modelSPIKE_conv, regressor_modelSPIKE, beta_out]=GenRegrModelTuk(dat,onset_SPIKE, SPIKE_template,tuksize)
%da Marta
if isempty(tuksize)
    win_tukey = 0.05;
else
    win_tukey=tuksize;
end;
win_tukey
regressor_modelSPIKE= zeros(1,length(dat));
%gh can be put OUT of loop since it remains fix
tukey = tukeywin(size(SPIKE_template,2), win_tukey)';   % Ho un 5% dei punti circa = 0, il resto � a 1
%  X =  [ detrend(SPIKE_template,1)'.*tukey'     ones(length(SPIKE_template),1)   detrend(1:length(SPIKE_template),0)'];  % detrend([1:length(ARTtemplate)].^2,0)' ];   %detrend(sin(1:length(ARTtemplate)),0)'];
X =  [ (SPIKE_template)'.*tukey'     zeros(length(SPIKE_template),1)  ];  % NO detrend and hence NO modelling of linear trends in template!!;
pinvX=pinv(X);
for rb = 1 : length(onset_SPIKE) %-  1                        % sottrazione ballisto per ballisto, di una media mobili di 100 ballisti successivi
    yy = [];
      yy = dat(onset_SPIKE(rb):onset_SPIKE(rb)+length(SPIKE_template)-1);
    beta(:,rb) = pinvX * yy;
    regressor_modelSPIKE(onset_SPIKE(rb):onset_SPIKE(rb)+length(SPIKE_template)-1) = X(:,1)*beta(1,rb);  %  se faccio X*beta(:,rb) tolgo proprio al media del pezzo e non va bene.
end

%-----------------------------------------------------
% CREIAMO UN REGRESSORE "MODELLO" CONVOLUTIVO
%-----------------------------------------------------

stick = zeros(1,length(dat));
stick(onset_SPIKE) = beta(1,:);
regressor_modelSPIKE_conv = conv(stick, SPIKE_template.*tukey);    %model2 = conv(stick, ARTtemplate_final{1}(1:256));
regressor_modelSPIKE_conv =   regressor_modelSPIKE_conv(1:min(  length(dat), length(regressor_modelSPIKE_conv) ));
beta_out=beta(1,:);
return;
