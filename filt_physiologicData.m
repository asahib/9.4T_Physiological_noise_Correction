function [pulse02_filt, pulse03_filt, resp02_filt, resp03_filt]=...
    filt_physiologicData(pulse02, pulse03, resp02, resp03,Hd)

pulse02_filt=filtfilt(Hd.Numerator,1,pulse02);
pulse03_filt=filtfilt(Hd.Numerator,1,pulse03);

resp02_filt=filtfilt(Hd.Numerator,1,resp02);
resp03_filt=filtfilt(Hd.Numerator,1,resp03);
end
