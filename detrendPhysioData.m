function detrendedPhysioData = detrendPhysioData(PhysioData, mintab)

yi=mintab(:,2);
xi=0:1:length(PhysioData)-1;
x=mintab(:,1);
y=interp1(x,yi,xi,'spline');
detrendedPhysioData=PhysioData(:)-y(:);

end