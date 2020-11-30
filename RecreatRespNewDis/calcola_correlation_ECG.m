function [ons_B, corART, M, corr_used, corr_vect, disp, num_peaks] = calcola_correlation_spike(dat, nsec, tot_nsec, ARTtemp, thresh, fc)

% - Calcola la correlazione.
% - Per ogni thresh_corr:
%   - Prima soglia su correlazione (thresh_corr) e su massimo (thresh.ampli) nel tratto lungo length(ARTtemp); indici -> firstpass
%   - Poi si selezionano per ogni picco due punti al lato del picco (distanza
%   tra indici > thresh.grad (gradient = met� (IMP) della distanza tra indici); indici -> Bons
%   - Si trova il max nell'intervallo tra i due indici -> ons_B


tART    = length(ARTtemp);
% se pi� di 100s velocizziamo spezzando in pezzi (bunches) da 100s.
if nsec > tot_nsec
    sub_nsec = [1:tot_nsec:nsec nsec];
else
    sub_nsec = [1 nsec];
end
tic; corART = [];  M = [];  
fprintf(['Calculating correlation. Bunch of seconds N. (over '  num2str(length(sub_nsec)-1)  '): '] )
for sub = 1 : length(sub_nsec)-1  % cos� se � � scalare viene 1 
    fprintf('\n')
    fprintf(['Bunch N. ' num2str(sub) '\n'])
    for t = 1 + fc*(sub_nsec(sub)-1) :  min([fc*(sub_nsec(sub+1)-1) length(dat)-tART])       %for t = 1 :  (fc)*nsec% length(dat)-length(ARTtemp)
        c = corrcoef(dat(t:t+tART-1),ARTtemp);
        corARTs(t-fc*(sub_nsec(sub)-1))=c(1,2);                                  % corART(t)=c(1,2);
        Ms(t-fc*(sub_nsec(sub)-1)) = max(abs(dat(t:t+round(tART*0.7)-1)));          % NEW round(tART/2), OLD tART              % M(t) = max(abs(dat(t:t+tART-1)));
        if    any(t == [1+fc*(sub_nsec(sub)-1): fc: fc*(sub_nsec(sub+1)+1)])                   %  any(t == [1: fc: fc*(nsec+1)])
            fprintf(['.'])
        end
    end
    corART = [corART corARTs];
    M = [M Ms];
    clear Ms corARTs
end;  toc

corr_vect = 0.5:0.01:1;   %= 0:0.01:1;
disp = zeros(1,length(corr_vect));
c = 0
for corr = corr_vect
    c = c + 1;
    thresh.corr = corr;
    clear firstpass Bons ons_BB
    %    firstpass=find(corART >  thresh.corr & M > thresh.ampli & M < 250 );%    % gh: Brutto cut-off per ecg e resp l'ho tolto
    firstpass=find(corART >  thresh.corr & M > thresh.ampli );
    
    ons_Ball{c} = [];
    if ~isempty(firstpass)
        Bons=[firstpass(1) firstpass(find(gradient(firstpass) > thresh.grad))    firstpass(end)   ];  % il gradiente speza in due la distanza  % ho 2 punti i lati del picco di correl. per ogni artefatto
        if mod(length(Bons),2)   % se ho un numero dispari di Bons
            fprintf('Numero dispari... \n')    % su claudia1 tresh_corr =0.45
            %figure, plot(diff(Bons))  % non si vede granch�
            %   keyboard             % Togliere o il primo o l'ultimo....
            disp(c) = 1 ;               % disp = 1 == problem
        end
        count = 0;
        for Q = 1:2:(length(Bons) - rem(length(Bons),2)  )  % -rem cos� funziona anche se dispari
            [Mas Indx] = max( corART(Bons(Q):Bons(Q+1)) );
            count = count + 1;
            ons_BB(count) = Bons(Q) + Indx - 1;
            clear Mas Indx
        end
        ons_Ball{c} = ons_BB;
    end;
    num_peaks(c) = length(ons_Ball{c});
end
%figure, subplot(3,1,1), plot(corr_vect, num_peaks,'*')
%hold on, plot(corr_vect,disp*5,'r')
%keyboard

% Cerchiamo la thresh_corr con cui otteniamo il max numero di picchi.
ind_nodisp = []; clear num_peaks_temp
fprintf('Searching for optimal correlation threshold...\n')
%tic;    
num_peaks_temp = num_peaks; 
while isempty(ind_nodisp)
    index_corr_max = find( num_peaks_temp == max( num_peaks_temp));  %possono essere anche pi� di uno; in tal caso... prendiamo il pi� piccolo
    ind_nodisp = find(disp(index_corr_max) == 0);               % cerchiamo gli indici in cui non abbiamo problemi (problema = num. dispari)
    num_peaks_temp(index_corr_max) = 0;     % azzeriamo i massimi, cos� cerchiamo nel loop succesivo i nuovi massimi
   fprintf('.')
   %t = toc;
  % if t>100
   %    fprintf('Not converging...')
    %   break
    %end
end
index_corr_max_nodisp = index_corr_max(ind_nodisp);   % controlliamo che per queste correl. non otteniamo 

corr_used = corr_vect(max(index_corr_max_nodisp))
ons_B =  ons_Ball{max(index_corr_max_nodisp)};  %...prendiamo il pi� grande

%figure,
%subplot(3,1,2),plot(corART)  %
%hold on, plot(ons_B,corART(ons_B),'r*')   %  plot(Bons,corART(Bons),'r*') 
%hold on, plot(corr_used*ones(1,length(corART)),'g')
%subplot(3,1,3), plot(dat)
%hold on, plot(ons_B,corART(ons_B),'r*')
%hold on, plot(thresh.ampli*ones(1,length(dat)),'g')

%fprintf('Check whether optimized, then press return... \n') 


