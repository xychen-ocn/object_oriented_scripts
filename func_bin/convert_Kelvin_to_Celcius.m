
%%%% function #2: convert Kelvin to Celcius:
function ds = convert_Kelvin_to_Celcius(ds, sstname)

lev1_fieldn = fieldnames(ds);
if ismember(sstname,lev1_fieldn)
    ds.(sstname) = ds.(sstname) - 273.15;
else
    if isstruct(ds.(lev1_fieldn{1}))
        lev2_fieldn = fieldnames(ds.(lev1_fieldn{1}));
        if ismember(sstname,lev2_fieldn)
            for i=1:length(lev1_fieldn)
                FN = lev1_fieldn{i};
                ds.(FN).(sstname) = ds.(FN).(sstname)-273.15;
            end
        end
    else
        disp('error!')
        return
    end
end

end

