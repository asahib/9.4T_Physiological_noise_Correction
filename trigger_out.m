function PhysioData_trigger_out = trigger_out(PhysioData)

data_count=0;

for i=1:length(PhysioData)-1
    if PhysioData(i)<4500
        data_count=data_count+1;
        PhysioData_trigger_out(data_count)=PhysioData(i);
    end 
end


end
