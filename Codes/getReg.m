function [MU,SIG] = getReg(data,param,setting,TT)

NN = size(TT,1);
[N,P] = size(data.Y);

ID = reshape(data.ID,[N*P,1]);

M = sum(ID==0);

T0 = data.T;

RR = zeros(N*P,1);
for p = 1:P
    RR((p-1)*N+1:p*N,1) = param.LAMBDA(p)^2*data.R(:,p);
end

if strcmp(setting.kernel,'SE')
    K_00 = exp(-exp(param.GAMMA)*(T0-T0').^2);
    K_T0 = exp(-exp(param.GAMMA)*(TT-T0').^2);
elseif strcmp(setting.kernel,'M15')
    K_00 = (1+sqrt(3).*exp(param.GAMMA).*abs(T0-T0')).*exp(-sqrt(3)*exp(param.GAMMA).*abs(T0-T0'));
    K_T0 = (1+sqrt(3).*exp(param.GAMMA).*abs(TT-T0')).*exp(-sqrt(3)*exp(param.GAMMA).*abs(TT-T0'));
elseif strcmp(setting.kernel,'M25')
    K_00 = (1+sqrt(5).*exp(param.GAMMA).*abs(T0-T0')+5/3*exp(2*param.GAMMA).*(T0-T0').^2).*exp(-sqrt(5)*exp(param.GAMMA).*abs(T0-T0'));
    K_T0 = (1+sqrt(5).*exp(param.GAMMA).*abs(TT-T0')+5/3*exp(2*param.GAMMA).*(TT-T0').^2).*exp(-sqrt(5)*exp(param.GAMMA).*abs(TT-T0'));
elseif strcmp(setting.kernel,'OU')
    K_00 = exp(-exp(param.GAMMA)*abs(T0-T0'));
    K_T0 = exp(-exp(param.GAMMA)*abs(TT-T0'));
end
KK_00 = kron(param.SIGMA,K_00) + diag(RR);

KK_00 = KK_00(ID==0,ID==0);


INV = KK_00\eye(M);


MU = zeros(NN,P);
NU = zeros(NN,P);

for n = 1:NN
    DD = kron(param.SIGMA,K_T0(n,:));
    
    DD = DD(:,ID==0);
    
    GG = reshape(data.MU-param.BIAS'-data.MEAN,[N*P,1]);
    GG = GG(ID==0);
    
    AA = DD*(INV*GG);
    MU(n,:) = AA' + param.BIAS';
    
    BB = param.SIGMA - DD*INV*DD';
    NU(n,:) = diag(BB)';
    
    GG = reshape(data.SIG,[N*P,1]);
    GG = GG(ID==0);
    
    CC = DD*INV*diag(GG);
    CC = CC*CC';
    NU(n,:) = NU(n,:) + diag(CC)';
end

SIG = sqrt(NU);


end