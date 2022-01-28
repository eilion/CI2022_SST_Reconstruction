function [data,Samples] = getSamples_Final(data,param,Samples,setting)

M = setting.nSamples;
P = length(data);

for p = 1:P
    TT = data(p).T0;
    
    RR = param.LAMBDA^2*ones(size(TT,1),1);
    
    AGE = data(p).T;
    N = size(AGE,1);
    
    if strcmp(setting.kernel,'SE')
        KK = exp(-exp(param.GAMMA(1))*(AGE-AGE').^2);
        KK_T0 = exp(-exp(param.GAMMA(1))*(TT-AGE').^2);
    elseif strcmp(setting.kernel,'M15')
        KK = (1+sqrt(3).*exp(param.GAMMA(1)).*abs(AGE-AGE')).*exp(-sqrt(3)*exp(param.GAMMA(1)).*abs(AGE-AGE'));
        KK_T0 = (1+sqrt(3).*exp(param.GAMMA(1)).*abs(TT-AGE')).*exp(-sqrt(3)*exp(param.GAMMA(1)).*abs(TT-AGE'));
    elseif strcmp(setting.kernel,'M25')
        KK = (1+sqrt(5).*exp(param.GAMMA(1)).*abs(AGE-AGE')+5/3*exp(2*param.GAMMA(1)).*(AGE-AGE').^2).*exp(-sqrt(5)*exp(param.GAMMA(1)).*abs(AGE-AGE'));
        KK_T0 = (1+sqrt(5).*exp(param.GAMMA(1)).*abs(TT-AGE')+5/3*exp(2*param.GAMMA(1)).*(TT-AGE').^2).*exp(-sqrt(5)*exp(param.GAMMA(1)).*abs(TT-AGE'));
    elseif strcmp(setting.kernel,'OU')
        KK = exp(-exp(param.GAMMA(1))*abs(AGE-AGE'));
        KK_T0 = exp(-exp(param.GAMMA(1))*abs(TT-AGE'));
    elseif strcmp(setting.kernel,'SE[SE]')
        GG = exp(-exp(param.GAMMA(1))*(AGE-AGE').^2);
        KK = (1+exp(param.GAMMA(2)).*(1-GG)).^(-1/2);
        
        GG_T0 = exp(-exp(param.GAMMA(1))*(TT-AGE').^2);
        KK_T0 = (1+exp(param.GAMMA(2)).*(1-GG_T0)).^(-1/2);
    end
    KK = KK + param.LAMBDA^2*diag(data(p).R);
    KK = data(p).SCALE^2.*KK;
    INV = KK\eye(N);
    
    KK_T0 = data(p).SCALE^2.*KK_T0;
    
    SIG = sqrt(data(p).SCALE^2*(1+RR)-sum(KK_T0.*(KK\KK_T0')',2));
    SIG = SIG.*ones(1,M);
    
    ID = ceil(size(Samples(p).MCMC_Samples,2)*rand(M,1));
    
    MU = data(p).BIAS + data(p).MEAN0 + KK_T0*(INV*(Samples(p).MCMC_Samples(:,ID)-data(p).BIAS-data(p).MEAN));
    
    Samples(p).SAMPLES_QUERY = normrnd(MU,SIG);
    
    data(p).LB0 = quantile(Samples(p).SAMPLES_QUERY,0.025,2);
    data(p).MD0 = quantile(Samples(p).SAMPLES_QUERY,0.5,2);
    data(p).UB0 = quantile(Samples(p).SAMPLES_QUERY,0.975,2);
end


end