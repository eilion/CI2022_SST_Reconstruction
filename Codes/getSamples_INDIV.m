function [data,Samples] = getSamples_INDIV(data,Samples,Model)

P = length(data);
nSamples = 10000;

for p = 1:P
    N = size(data(p).Y,1);
    
    Samples(p).MCMC_Samples_INDIV = zeros(N,nSamples);
end

% Get individual samples:
parfor p = 1:P
    Samples(p).MCMC_Samples_INDIV = MCMC_Sampling_INDIV(data(p),Model,nSamples);
    
    data(p).LB_INDIV = quantile(Samples(p).MCMC_Samples_INDIV,0.025,2);
    data(p).MD_INDIV = quantile(Samples(p).MCMC_Samples_INDIV,0.5,2);
    data(p).UB_INDIV = quantile(Samples(p).MCMC_Samples_INDIV,0.975,2);
end


end