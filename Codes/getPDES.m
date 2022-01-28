function [PDEVS] = getPDES(data,param,TT,YY,setting)

PDEVS = struct('LOGLIK',cell(1,1));

PDEVS.LOGLIK = 0;
PDEVS.SCALE = 0;
PDEVS.BIAS = 0;
PDEVS.GAMMA = zeros(2,1);
PDEVS.LAMBDA = 0;


N = size(TT,1);

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

ALPHA = KK_00\(YY-data.BIAS-data.MEAN);

PDEVS.LOGLIK = PDEVS.LOGLIK - 0.5.*data.SCALE.^(-2).*(YY-data.BIAS-data.MEAN)'*ALPHA - sum(log(diag(chol(KK_00)))) - N*log(data.SCALE);

% SCALE and BIAS:
PDEVS.SCALE = PDEVS.SCALE + data.SCALE.^(-3).*(YY-data.BIAS-data.MEAN)'*ALPHA - N*data.SCALE.^(-1);
PDEVS.BIAS = PDEVS.BIAS + data.SCALE.^(-2).*sum(ALPHA);

QQ = (data.SCALE.^(-2)*(ALPHA*ALPHA')-KK_00\eye(N))';

% GAMMA:
if strcmp(setting.kernel,'SE[SE]')
    PK = - 0.5.*K_00.^3.*exp(param.GAMMA(2)).*G_00.*exp(param.GAMMA(1)).*(TT-TT').^2;
    PDEVS.GAMMA(1) = PDEVS.GAMMA(1) + 0.5*sum(sum(QQ.*PK));
    
    % PK = - 0.5.*K_00.^3.*exp(param.GAMMA(2)).*(1-G_00);
    % PK = PK.*data.SCALE.^(-2);
    % PDEVS.GAMMA(2) = PDEVS.GAMMA(2) + 0.5*sum(sum(QQ.*PK));
else
    if strcmp(setting.kernel,'SE')
        PK = K_00.*(-exp(param.GAMMA(1))*(TT-TT').^2);
    elseif strcmp(setting.kernel,'M15')
        PK = - 3.*exp(param.GAMMA(1)).*abs(TT-TT').^2.*exp(-sqrt(3)*exp(param.GAMMA(1)).*abs(TT-TT'));
        PK = PK.*exp(param.GAMMA(1));
    elseif strcmp(setting.kernel,'M25')
        PK = - (5/3.*exp(param.GAMMA(1)).*abs(TT-TT').^2+5*sqrt(5)/3.*exp(2*param.GAMMA(1)).*abs(TT-TT').^3).*exp(-sqrt(5)*exp(param.GAMMA(1)).*abs(TT-TT'));
        PK = PK.*exp(param.GAMMA(1));
    elseif strcmp(setting.kernel,'OU')
        PK = K_00.*(-exp(param.GAMMA(1))*abs(TT-TT'));
    end
    PDEVS.GAMMA(1) = PDEVS.GAMMA(1) + 0.5*sum(sum(QQ.*PK));
end

% LAMBDA:
PK = 2*param.LAMBDA*diag(data.R);
PDEVS.LAMBDA = PDEVS.LAMBDA + 0.5*sum(sum(QQ.*PK));


end