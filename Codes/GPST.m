function [data,param,timeline] = GPST(data,param,Model,setting)

beta1 = 0.9;
beta2 = 0.999;
epsilon = 1e-8;
gamma = 1e-3;

max_iters = setting.max_iters;

[N,P] = size(data.Y);


timeline = struct('name',cell(P,1),'TT',cell(P,1),'MU',cell(P,1),'SIG',cell(P,1));
for p = 1:P
    timeline(p).name = ['signal ',num2str(p)];
    
    timeline(p).TT = linspace(setting.st,setting.ed,setting.NT)';
    timeline(p).MEAN = interp1(data.MEAN0(:,1),data.MEAN0(:,2),timeline(p).TT);
    
    timeline(p).MU = zeros(setting.NT,floor(max_iters/1000));
    timeline(p).SIG = zeros(setting.NT,floor(max_iters/1000));
end


Mw_SIGMA = zeros(P);
Vw_SIGMA = zeros(P);

Mw_GAMMA = 0;
Vw_GAMMA = 0;

Mw_LAMBDA = zeros(P,1);
Vw_LAMBDA = zeros(P,1);

Mw_BIAS = zeros(P,1);
Vw_BIAS = zeros(P,1);

Mw_MU = zeros(N,P);
Vw_MU = zeros(N,P);

Mw_SIG = zeros(N,P);
Vw_SIG = zeros(N,P);


LOGLIK = zeros(max_iters,1);


for r = 1:max_iters
    
    %{
    if strcmp(setting.var,'heteroscedastic')
        if rem(r,2000) == 1 && r > 1
            [data,param] = updateVar(data,param,setting);
        end
    end
    %}
    
    
    PDEV_SIGMA = zeros(P,P);
    PDEV_GAMMA = 0;
    PDEV_LAMBDA = zeros(P,1);
    
    PDEV_BIAS = zeros(P,1);
    
    PDEV_MU = zeros(N,P);
    PDEV_SIG = zeros(N,P);
    
    EPS_X = normrnd(0,1,[N,P]);
    
    
    % p(Y|X):
    for p = 1:P
        PDEVS = getPDES_Y(data,Model,EPS_X(:,p),p);
        
        PDEV_MU(:,p) = PDEV_MU(:,p) + PDEVS.MU;
        PDEV_SIG(:,p) = PDEV_SIG(:,p) + PDEVS.SIG;
        
        LOGLIK(r) = LOGLIK(r) + PDEVS.LOGLIK;
    end
    

    % p(X):
    PDEVS = getPDES(data,param,EPS_X,setting);
    
    PDEV_SIGMA = PDEV_SIGMA + PDEVS.SIGMA;
    PDEV_GAMMA = PDEV_GAMMA + PDEVS.GAMMA;
    PDEV_LAMBDA = PDEV_LAMBDA + PDEVS.LAMBDA;
    
    PDEV_BIAS = PDEV_BIAS + PDEVS.BIAS;
    
    PDEV_MU = PDEV_MU + PDEVS.MU;
    PDEV_SIG = PDEV_SIG + PDEVS.SIG;
    
    LOGLIK(r) = LOGLIK(r) + PDEVS.LOGLIK;
    
    
    % q(X,X0):
    PDEV_SIG = PDEV_SIG + 1./data.SIG;
    PDEV_SIG(isnan(PDEV_SIG)) = 0;
    
    LOGLIK(r) = LOGLIK(r) + sum(log(data.SIG(data.ID==0)));
    
    
    Mw_BIAS = beta1*Mw_BIAS - (1-beta1)*PDEV_BIAS;
    Vw_BIAS = beta2*Vw_BIAS + (1-beta2)*PDEV_BIAS.*PDEV_BIAS;
    
    Mw_SIGMA = beta1*Mw_SIGMA - (1-beta1)*PDEV_SIGMA;
    Vw_SIGMA = beta2*Vw_SIGMA + (1-beta2)*PDEV_SIGMA.*PDEV_SIGMA;
    
    Mw_GAMMA = beta1*Mw_GAMMA - (1-beta1)*PDEV_GAMMA;
    Vw_GAMMA = beta2*Vw_GAMMA + (1-beta2)*PDEV_GAMMA.*PDEV_GAMMA;
    
    Mw_LAMBDA = beta1*Mw_LAMBDA - (1-beta1)*PDEV_LAMBDA;
    Vw_LAMBDA = beta2*Vw_LAMBDA + (1-beta2)*PDEV_LAMBDA.*PDEV_LAMBDA;
    
    Mw_MU = beta1*Mw_MU - (1-beta1)*PDEV_MU;
    Vw_MU = beta2*Vw_MU + (1-beta2)*PDEV_MU.*PDEV_MU;
    
    Mw_SIG = beta1*Mw_SIG - (1-beta1)*PDEV_SIG;
    Vw_SIG = beta2*Vw_SIG + (1-beta2)*PDEV_SIG.*PDEV_SIG;
    
    
    % param.SIGMA_L = param.SIGMA_L - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw_SIGMA./(sqrt(Vw_SIGMA)+epsilon);
    % param.SIGMA = param.SIGMA_L*param.SIGMA_L' + 1e-6*eye(P);
    
    % param.GAMMA = param.GAMMA - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw_GAMMA./(sqrt(Vw_GAMMA)+epsilon);
    % param.LAMBDA = param.LAMBDA - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw_LAMBDA./(sqrt(Vw_LAMBDA)+epsilon);
    
    % param.BIAS = param.BIAS - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw_BIAS./(sqrt(Vw_BIAS)+epsilon);
    
    data.MU = data.MU - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw_MU./(sqrt(Vw_MU)+epsilon);
    data.SIG = data.SIG - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw_SIG./(sqrt(Vw_SIG)+epsilon);
    
    
    
    if rem(r,100) == 0
        disp(['   Iterations ',num2str(r),'/',num2str(max_iters),' are done.']);
    end
    
    
    if rem(r,1000) == 0
        
        TT = (timeline(p).TT-setting.mean_T)/setting.stdv_T;
        
        [MU,SIG] = getReg(data,param,setting,TT);
        
        if strcmp(setting.var,'heteroscedastic')
            RR = getVar(data,param,setting,TT);
        elseif strcmp(setting.var,'homoscedastic')
            RR = zeros(size(MU,1),P);
            for p = 1:P
                RR(:,p) = param.LAMBDA(p)^2*mean(data.R(data.ID(:,p)==0,p));
            end
        end
        
        for p = 1:P
            timeline(p).MU(:,r/1000) = MU(:,p) + timeline(p).MEAN;
            timeline(p).SIG(:,r/1000) = sqrt(SIG(:,p).^2+RR(:,p));
        end
        

        close all;
        
        fig = figure;
        
        for p = 1:P
            subplot(P,1,p);
            hold on;
            
            UB = timeline(p).MU(:,r/1000) + 1.96*timeline(p).SIG(:,r/1000);
            LB = timeline(p).MU(:,r/1000) - 1.96*timeline(p).SIG(:,r/1000);
            MB = timeline(p).MU(:,r/1000);
            
            T0 = data.T*setting.stdv_T + setting.mean_T;
            % Y0 = (data.Y(:,p)-0.044)/0.033;
            
            for n = 1:size(T0,1)
                plot([T0(n),T0(n)],[data.LB(n,p),data.UB(n,p)],'g','LineWidth',2);
                % plot([T0(n),T0(n)],[Y0(n)-2.5,Y0(n)+2.5],'g','LineWidth',2);
            end
            
            plot(timeline(p).TT,LB,'k','LineWidth',2);
            plot(timeline(p).TT,MB,'--k','LineWidth',2);
            plot(timeline(p).TT,UB,'k','LineWidth',2);
            
            ID = isnan(data.MU(:,p));
            plot(T0(~ID),data.MU(~ID,p)-1.96*data.SIG(~ID,p),'r','LineWidth',2);
            plot(T0(~ID),data.MU(~ID,p)+1.96*data.SIG(~ID,p),'r','LineWidth',2);
            
            % plot(T0(~ID),Y0(~ID),'*b','LineWidth',2);
            plot(T0(~ID),data.MD(~ID,p),'*b','LineWidth',2);
            
            xlim([timeline(p).TT(1) timeline(p).TT(end)]);
            ylim([-3 32]);
        end
        
        set(fig,'Position',[20 20 1200 50*(P-1)+350*P]);
        movegui(fig,'center');
        drawnow;
    end
end


end