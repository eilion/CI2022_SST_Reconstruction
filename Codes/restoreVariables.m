function [data] = restoreVariables(data,setting)

P = length(data);

for p = 1:P
    data(p).T = data(p).T*setting.stdv_T + setting.mean_T;
    data(p).T0 = data(p).T0*setting.stdv_T + setting.mean_T;
end


end