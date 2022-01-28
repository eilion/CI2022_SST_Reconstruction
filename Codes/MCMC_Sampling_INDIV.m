function [SAMPLES] = MCMC_Sampling_INDIV(data,Model,nSamples)

UK37 = data.Y;

N = size(UK37,1);

SAMPLES = NaN*ones(N,nSamples);

BIAS = (UK37-0.044)/0.033;

SAM = BIAS;

MM = interp1(Model(:,1),Model(:,2),SAM);
SS = interp1(Model(:,1),Model(:,3),SAM);

MIN = Model(3,1);
MAX = Model(end-2,1);

LOGLIK_old = - 0.5*((UK37-MM)./SS).^2 - log(SS) - 0.5*log(2*pi);

MAX_ITERS = 1000 + 100*(nSamples-1);

n = 0;
for r = 1:MAX_ITERS
    if rem(r,100) == 1
        rand_seed = log(rand(N,100));
        t = 1;
    end
    
    SAM_new = SAM + normrnd(0,1,[N,1]);
    
    MM = interp1(Model(:,1),Model(:,2),SAM_new);
    SS = interp1(Model(:,1),Model(:,3),SAM_new);
    
    LOGLIK_new = - 0.5*((UK37-MM)./SS).^2 - log(SS) - 0.5*log(2*pi);
    LOGLIK_new((SAM_new<MIN)|(SAM_new>MAX)) = -inf;
    
    ID = (LOGLIK_new-LOGLIK_old>rand_seed(:,t));
    
    SAM(ID) = SAM_new(ID);
    LOGLIK_old(ID) = LOGLIK_new(ID);
    
    t = t + 1;
    
    if r == 1000
        n = n + 1;
        SAMPLES(:,n) = SAM;
    end
    
    if r > 1000 && rem(r-1000,100) == 0
        n = n + 1;
        SAMPLES(:,n) = SAM;
    end
end


end