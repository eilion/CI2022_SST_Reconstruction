function [data,Samples] = getSamples(data,param,Samples,Model,setting)

P = length(data);
nSamples = 1000;

for p = 1:P
    N = size(data(p).Y,1);
    
    Samples(p).MCMC_Samples = zeros(N,nSamples);
end

% Get individual samples:
parfor p = 1:P
    Samples(p).MCMC_Samples = MCMC_Sampling(data(p),param,Model,nSamples,setting);
    
    data(p).LB = quantile(Samples(p).MCMC_Samples,0.025,2);
    data(p).MD = quantile(Samples(p).MCMC_Samples,0.5,2);
    data(p).UB = quantile(Samples(p).MCMC_Samples,0.975,2);
end


end