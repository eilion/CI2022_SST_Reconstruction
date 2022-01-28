function [SAMPLES] = MCMC_Sampling(data,param,Model,nSamples,setting)

epsilon = 1e-1;
L = 10;

UK37 = data.Y;
TT = data.T;

N = size(UK37,1);

if strcmp(setting.kernel,'SE')
    K_00 = exp(-exp(param.GAMMA(1))*(TT-TT').^2);
elseif strcmp(setting.kernel,'M15')
    K_00 = (1+sqrt(3).*exp(param.GAMMA(1)).*abs(TT-TT')).*exp(-sqrt(3)*exp(param.GAMMA(1)).*abs(TT-TT'));
elseif strcmp(setting.kernel,'M25')
    K_00 = (1+sqrt(5).*exp(param.GAMMA(1)).*abs(TT-TT')+5/3*exp(2*param.GAMMA(1)).*(TT-TT').^2).*exp(-sqrt(5)*exp(param.GAMMA(1)).*abs(TT-TT'));
elseif strcmp(setting.kernel,'OU')
    K_00 = exp(-exp(param.GAMMA(1))*abs(TT-TT'));
elseif strcmp(setting.kernel,'SE[SE]')
    G_00 = exp(-exp(param.GAMMA(1))*(TT-TT').^2);
    K_00 = (1+exp(param.GAMMA(2)).*(1-G_00)).^(-1/2);
end
KK_00 = K_00 + param.LAMBDA^2*diag(data.R);
KK_00 = data.SCALE^2.*KK_00;

INV = KK_00\eye(N);


SAMPLES = NaN*ones(N,nSamples);


BIAS = mean((UK37-0.044)/0.033);


SAM = BIAS*ones(N,1);
PHI = normrnd(0,1,[N,1]);

MM = interp1(Model(:,1),Model(:,2),SAM);
SS = interp1(Model(:,1),Model(:,3),SAM);

LOGLIK_old = sum(-0.5*PHI.^2-0.5*log(2*pi));
LOGLIK_old = LOGLIK_old + sum(-0.5*((UK37-MM)./SS).^2-log(SS)-0.5*log(2*pi));
LOGLIK_old = LOGLIK_old - 0.5*(SAM-data.BIAS-data.MEAN)'*INV*(SAM-data.BIAS-data.MEAN);


MAX_ITERS = 10000 + 200*(nSamples-1);

rand_seed = log(rand(MAX_ITERS,1));

for r = 1:MAX_ITERS
    
    MM = interp1(Model(:,1),Model(:,2),SAM);
    SS = interp1(Model(:,1),Model(:,3),SAM);
    
    PM = interp1(Model(:,1),Model(:,4),SAM);
    PS = interp1(Model(:,1),Model(:,5),SAM);
    
    PDEV = - INV*(SAM-data.BIAS-data.MEAN) + (UK37-MM).*SS.^(-2).*PM + (UK37-MM).^2.*SS.^(-3).*PS - PS./SS;
    
    PHI_new = PHI + 0.5*epsilon*PDEV;
    SAM_new = SAM;
    
    for ll = 1:L-1
        SAM_new = SAM_new + epsilon*PHI_new;
        
        MM = interp1(Model(:,1),Model(:,2),SAM_new);
        SS = interp1(Model(:,1),Model(:,3),SAM_new);
        
        PM = interp1(Model(:,1),Model(:,4),SAM_new);
        PS = interp1(Model(:,1),Model(:,5),SAM_new);
        
        PDEV = - INV*(SAM_new-data.BIAS-data.MEAN) + (UK37-MM).*SS.^(-2).*PM + (UK37-MM).^2.*SS.^(-3).*PS - PS./SS;
        PHI_new = PHI_new + epsilon*PDEV;
    end
    
    SAM_new = SAM_new + epsilon*PHI_new;
    
    MM = interp1(Model(:,1),Model(:,2),SAM_new);
    SS = interp1(Model(:,1),Model(:,3),SAM_new);
    
    PM = interp1(Model(:,1),Model(:,4),SAM_new);
    PS = interp1(Model(:,1),Model(:,5),SAM_new);
    
    PDEV = - INV*(SAM_new-data.BIAS-data.MEAN) + (UK37-MM).*SS.^(-2).*PM + (UK37-MM).^2.*SS.^(-3).*PS - PS./SS;
    PHI_new = PHI_new + 0.5*epsilon*PDEV;
    
    
    MM = interp1(Model(:,1),Model(:,2),SAM_new);
    SS = interp1(Model(:,1),Model(:,3),SAM_new);
    
    LOGLIK_new = sum(-0.5*PHI_new.^2-0.5*log(2*pi));
    LOGLIK_new = LOGLIK_new + sum(-0.5*((UK37-MM)./SS).^2-log(SS)-0.5*log(2*pi));
    LOGLIK_new = LOGLIK_new - 0.5*(SAM_new-data.BIAS-data.MEAN)'*INV*(SAM_new-data.BIAS-data.MEAN);
    
    
    if LOGLIK_new - LOGLIK_old > rand_seed(r)
        SAM = SAM_new;
        PHI = PHI_new;
        LOGLIK_old = LOGLIK_new;
    end
    
    LOGLIK_old = LOGLIK_old - sum(-0.5*PHI.^2-0.5*log(2*pi));
    PHI = normrnd(0,1,[N,1]);
    LOGLIK_old = LOGLIK_old + sum(-0.5*PHI.^2-0.5*log(2*pi));
    
    if r == 10000
        t = 1;
        SAMPLES(:,t) = SAM;
    end
    
    if r > 10000 && rem(r-10000,200) == 0
        t = t + 1;
        SAMPLES(:,t) = SAM;
    end
end


end