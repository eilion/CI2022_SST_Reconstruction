function [data,param] = learnParameters(data,param,setting)

beta1 = 0.9;
beta2 = 0.999;
epsilon = 1e-8;
gamma = 1e-3;

max_iters = setting.max_iters;

P = length(data);

YY = cell(P,1);
TT = cell(P,1);

for p = 1:P
    TT{p} = data(p).T;
    YY{p} = data(p).MD_INDIV;
end


Mw_SCALE = zeros(P,1);
Vw_SCALE = zeros(P,1);

Mw_BIAS = zeros(P,1);
Vw_BIAS = zeros(P,1);

Mw_GAMMA = zeros(2,1);
Vw_GAMMA = zeros(2,1);

Mw_LAMBDA = 0;
Vw_LAMBDA = 0;

LOGLIK = zeros(max_iters,1);

for r = 1:max_iters
    
    PDEV_SCALE = zeros(P,1);
    PDEV_BIAS = zeros(P,1);
    
    PDEV_GAMMA = zeros(2,1);
    PDEV_LAMBDA = 0;
    
    LGLK = 0;
    
    parfor p = 1:P
        PDEVS = getPDES(data(p),param,TT{p},YY{p},setting);
        
        PDEV_SCALE(p) = PDEV_SCALE(p) + PDEVS.SCALE;
        PDEV_BIAS(p) = PDEV_BIAS(p) + PDEVS.BIAS;
        
        PDEV_GAMMA = PDEV_GAMMA + PDEVS.GAMMA;
        PDEV_LAMBDA = PDEV_LAMBDA + PDEVS.LAMBDA;
        
        LGLK = LGLK + PDEVS.LOGLIK;
    end
    
    LOGLIK(r) = LGLK;
    
    
    Mw_SCALE = beta1*Mw_SCALE - (1-beta1)*PDEV_SCALE;
    Vw_SCALE = beta2*Vw_SCALE + (1-beta2)*PDEV_SCALE.*PDEV_SCALE;
    
    Mw_BIAS = beta1*Mw_BIAS - (1-beta1)*PDEV_BIAS;
    Vw_BIAS = beta2*Vw_BIAS + (1-beta2)*PDEV_BIAS.*PDEV_BIAS;
    
    Mw_GAMMA(1) = beta1*Mw_GAMMA(1) - (1-beta1)*PDEV_GAMMA(1);
    Vw_GAMMA(1) = beta2*Vw_GAMMA(1) + (1-beta2)*PDEV_GAMMA(1).*PDEV_GAMMA(1);
    
    Mw_LAMBDA = beta1*Mw_LAMBDA - (1-beta1)*PDEV_LAMBDA;
    Vw_LAMBDA = beta2*Vw_LAMBDA + (1-beta2)*PDEV_LAMBDA.*PDEV_LAMBDA;
    
    for p = 1:P
        data(p).SCALE = abs(data(p).SCALE - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw_SCALE(p)./(sqrt(Vw_SCALE(p))+epsilon));
        data(p).BIAS = data(p).BIAS - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw_BIAS(p)./(sqrt(Vw_BIAS(p))+epsilon);
    end
    
    param.GAMMA = param.GAMMA - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw_GAMMA./(sqrt(Vw_GAMMA)+epsilon);
    param.LAMBDA = abs(param.LAMBDA - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw_LAMBDA./(sqrt(Vw_LAMBDA)+epsilon));
end

param.LOGLIK = LOGLIK;


end