function tseries=get_tseries(Vol,mask)

mask_tseries=reshape(Vol, size(Vol,1)*size(Vol,2)*size(Vol,3), size(Vol,4));
tseries=mask_tseries(mask>0,:);
end
