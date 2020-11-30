function [SPIKE_template, onset_SPIKE]=FindFirstTemplate(dat,thresh,Fc,varargin)

%please define varargin=[firstx secondx];


%da Marta  sistemato anche per respirazione
%variabili necessari:
%per ECG:
   % thresh.grad =  mean(diff(ntrig)*1/400);  % distanza fra trigger
   % thresh.ampli = 2048;  %  half the dynamic range
   % Fc=400;


   % chose = 'n';
   chose=nargin<4;
    if chose == 1
        figure, plot(dat)
        zoom xon
        tp=0;
        while tp==0 
            tp=waitforbuttonpress;
        end
        [xb yb] = ginput(2);
        zoom off
        close(gcf)
    else
        xb=varargin{1};
    end


    %xb = xbb{i};   
    firstARTtemp=dat( round(xb(1)):round(xb(2)) );  %12900    12996
%    figure(1), plot(firstARTtemp)
    Info_first.ARTtemplate = firstARTtemp;
    Info_first.xb = xb; 
    
  %-----------------------------------------------------
    % CALCOLIAMO LA CORRELAZIONE E TROVIAMO ALMENO 20 SPIKE SUL canale considerato (O2) 
    %-----------------------------------------------------
    fprintf('Looking for 20 spikes (at least) by correlation with a first template... ')
    clear ARTtemp tART nsec tot_nsec sub_nsec corARTs Ms corART M 
    ARTtemp = firstARTtemp;
    nsec = 80;
    tot_nsec = 40;    % se si vuole velocizzare dare tot_nsec < di nsec, altrimenti tot_nsec= nsec;
    clear ons_B corART M corr_used
    [ons_B, corART, M, corr_used, corr_vect, disp, num_peaks] = calcola_correlation_ECG(dat, nsec, tot_nsec, ARTtemp, thresh, Fc);
    Info_first20.thresh = thresh; 
    Info_first20.thresh.corr = corr_used;
    Info_first20.ons_B = ons_B;
    Info_first20.corART = corART;
    Info_first20.M = M;
    Info_first20.corr_vect = corr_vect;
    Info_first20.disp = disp;
    Info_first20.num_peaks = num_peaks;
    clear thresh
    close(gcf)
  
    
       %-----------------------------------------------------
    % CREIAMO UN TEMPLATE MEDIO 
    %gh molti param erano 'hard coded' quindi ho modificato assai
    %-----------------------------------------------------
    fprintf('Creating a preliminary template by averaging 20% of equally spaced spikes found...\n ')

   %gh find local correlation maximum if several correlation maxima are
   %found for same peak

   idx=find(diff(ons_B)>length(ARTtemp)/2);
   if (any(diff(idx)-1))
       [Y,newons_B(1)]=max(Info_first20.corART(1:ons_B(idx(1))));
       for k=2:length(idx)-1
           [Y,I]=max(Info_first20.corART(ons_B(idx(k-1)):ons_B(idx(k))));
           newons_B(k)=I+ons_B(idx(k-1));
       end;
       [Y,I]=max(Info_first20.corART(ons_B(idx(end)):end));
       newons_B(length(idx)+1)=I+ons_B(idx(end));
       ons_B=newons_B;
   end;
   ons_B=ons_B(ons_B>0);
   salta=ceil(length(ons_B)/20);
   tot_ART = ons_B(1:salta:end);   % gh lunghezza specificata in precedenza> length(ons_Bchan)
    ARTtemplate = [];     %cell(1,length(ons_Bchan));   % = []
    for in = 1  %length(ons_B)-20-1  % prendiamo gli ultimi 20 di quelli trovati (a volte i primi sembrano schifezze
        ARTtmp=[];
        for rb = tot_ART
            ARTtmp=[ARTtmp; dat(rb:rb+length(ARTtemp)-1)];   % tART(in)
        end
        ARTtemplate{in}=median(ARTtmp);         %  in-in+1      
    end 
    Info_first20.ARTtemplate = ARTtemplate{1};
%    figure(1), hold on,  plot(ARTtemplate{1},'r')
%    hold on, plot(firstARTtemp,'b')
    
    
    
    %-----------------------------------------------------
    % OPTIMIMIZING THE CORRELATION THRESHOLD
    %-----------------------------------------------------    


    fprintf('Looking for all spikes by correlation with the preliminary template...\n ')
    clear ARTtemp tART nsec tot_nsec sub_nsec corARTs Ms corART M 
    ARTtemp = ARTtemplate{1};%use new template
    nsec = floor(length(dat)/Fc)-3;    % -1 perch� altrim l'ultimo non riesce a farlo 
    tot_nsec = 25;    % se si vuole velocizzare dare tot_nsec < di nsec, altrimenti tot_nsec= nsec;
   % thresh.grad = round(85/4);    %  85 per 1024, /4 per 256      % 80 per il BKG channel;          200 per il ballisto sui canali;
   % thresh.ampli = 50;  %   100 per il BKG channel;          60 per il ballisto sui canali;
    thresh.grad=Info_first20.thresh.grad;
    thresh.ampli=min([Info_first20.thresh.ampli max(ARTtemp)/2]);
    clear ons_B corART M corr_used
    [ons_B, corART, M, corr_used, corr_vect, disp, num_peaks] = calcola_correlation_ECG(dat, nsec, tot_nsec, ARTtemp, thresh, Fc);
    Info_morerobust.thresh = thresh; 
    Info_morerobust.thresh.corr = corr_used;
    Info_morerobust.ons_B = ons_B;
    Info_morerobust.corART = corART;
    Info_morerobust.M = M;
    Info_morerobust.corr_vect = corr_vect;
    Info_morerobust.disp = disp;
    Info_morerobust.num_peaks = num_peaks;
    clear thresh
 %   close(gcf)
    %-----------------------------------------------------
    % CREIAMO UN TEMPLATE MEDIO PER tutti gli SPIKE trovati sul canale considerato 
    % lungo come il template esterno (150 punti)
    %-----------------------------------------------------
    
  % keyboard
   %gh find local correlation maximum if several correlation maxima are
   %found for same peak

    idx=find(diff(ons_B)>length(ARTtemp)/2);
    if (any(diff(idx)-1))
        clear newons_B;
        [Y,newons_B(1)]=max(Info_morerobust.corART(1:ons_B(idx(1))));
        for k=2:length(idx)-1
            [Y,I]=max(Info_morerobust.corART(ons_B(idx(k-1)):ons_B(idx(k))));
            newons_B(k)=I+ons_B(idx(k-1));
        end;
        [Y,I]=max(Info_morerobust.corART(ons_B(idx(end)):end));
        newons_B(length(idx)+1)=I+ons_B(idx(end));

        ons_B=newons_B;
    end;
    fprintf('Creating a more robust template by averaging all the spikes found...\n ')
    clear tot_ART  ARTtemplate_final
    ons_B=ons_B(ons_B>0);
    tot_ART = length(ons_B);   % > length(ons_Bchan)
    length_temp =  length(firstARTtemp);  % 2 s


    % delay = 300;  %>0 se voglio anticipare
    ARTtemplate = [];     %cell(1,length(ons_Bchan));   % = []
    for in = 1 %: length(ons_B) - tot_ART + 1                     % faccio media di ballisti non su tutto il tracciato ma ogni 100 ballisti
        ARTtmp=[];
        for rb = in : in + tot_ART -1       %  1 : length(ons_Bchan) %- 1     %did exceed length trial ...
            ARTtmp=[ARTtmp; dat(ons_B(rb):ons_B(rb)+length_temp-1)];   % tART(in)
        end
        ARTtemplate_final{in}=median(ARTtmp);                  %mean o median
    end    
    Info_morerobust.ARTtemplate = ARTtemplate_final{1};
    %Info_morerobust.delay = delay;
    %figure(1), hold on,  plot(ARTtemplate_final{1},'g')

    %init_spike_template = ARTtemplate_final{1};
    
      % ITERIAMO 1 ULTIMA VOLTA
    %-----------------------------------------------------
    % RI-CALCOLIAMO LA CORRELAZIONE E TROVIAMO TUTTI GLI ARTEFATTI SUL canale
    % considerato a partire da ARTtemplate_final
    % OPTIMIMIZING THE CORRELATION THRESHOLD
    %-----------------------------------------------------
    fprintf('Looking for all spikes by correlation with the more robust template...\n ')
    clear ARTtemp tART nsec tot_nsec sub_nsec corARTs Ms corART M 
    %ARTtemp = init_spike_template; %ARTtemplate_final{1};
    ARTtemp = ARTtemplate_final{1};
    nsec = ceil(length(dat)/Fc);    % gh to really get the last ones in trace... but may create future problems in correlation calc..
    tot_nsec = 25;    % se si vuole velocizzare dare tot_nsec < di nsec, altrimenti tot_nsec= nsec;
    %tot_nsec=nsec
    thresh.grad=Info_morerobust.thresh.grad;
    thresh.ampli=min([Info_morerobust.thresh.ampli max(ARTtemp)/2]);
    %thresh.grad = round(120/4);    %  85/4;  se i = 8 120/4
    %thresh.ampli = 50;  %   100 per il BKG channel;          60 per il ballisto sui canali;
    
    clear ons_B corART M corr_used
    [ons_B, corART, M, corr_used, corr_vect, disp, num_peaks] = calcola_correlation_ECG(dat, nsec, tot_nsec, ARTtemp, thresh, Fc);
    Info_final.thresh = thresh; 
    Info_final.thresh.corr = corr_used;
    Info_final.ons_B = ons_B;
    Info_final.corART = corART;
    Info_final.M = M;
    Info_final.corr_vect = corr_vect;
    Info_final.disp = disp;
    Info_final.num_peaks = num_peaks;
    clear thresh
    
    %-----------------------------------------------------
    % CREIAMO UN TEMPLATE MEDIO PER tutti gli SPIKE trovati sul canale considerato 
    %-----------------------------------------------------
    fprintf('Creating a final template by averaging all the spikes found... ')
    clear tot_ART  ARTtemplate_final_final
    
    %gh find local correlation maximum

    idx=find(diff(ons_B)>length(ARTtemp)/2);
    if (any(diff(idx)-1))
        clear newons_B ;
        [Y,newons_B(1)]=max(Info_final.corART(1:ons_B(idx(1))));
        for k=2:length(idx)
            [Y,I]=max(Info_final.corART(ons_B(idx(k-1)):ons_B(idx(k))));
            newons_B(k)=I+ons_B(idx(k-1));
            clear Y I;
        end;
        [Y,I]=max(Info_final.corART(ons_B(idx(end)):end));
        newons_B(length(idx)+1)=I+ons_B(idx(end));
        ons_B=newons_B;
    end;

    tot_ART = length(ons_B);   % > length(ons_Bchan)
    length_temp = length(ARTtemp);%L;  %fc*2; %  length(firstARTtemp);  % 2 s
    
    ARTtemplate = [];     %cell(1,length(ons_Bchan));   % = []
    for in = 1 %: length(ons_B) - tot_ART + 1                     % faccio media di ballisti non su tutto il tracciato ma ogni 100 ballisti
        ARTtmp=[];
        for rb = in : in + tot_ART -1       %  1 : length(ons_Bchan) %- 1     %did exceed length trial ...
            ARTtmp=[ARTtmp; dat(ons_B(rb):ons_B(rb)+length_temp-1)];   % tART(in)
        end
        ARTtemplate_final_final{in}=median(ARTtmp);                  %mean o median
    end    
    Info_final.ARTtemplate = ARTtemplate_final_final{1};
    
  %  figure(1), hold on,  plot(ARTtemplate_final_final{1},'m')
    
    %figure, plot(ARTtemplate_final_final{1})
    %hold on, plot(firstARTtemp,'r')
    %plot(ARTtemplate_final_final{1},'g')
    
    SPIKE_template = ARTtemplate_final_final{1};
    onset_SPIKE = ons_B;
  

return;

