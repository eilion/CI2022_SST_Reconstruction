function [data,param,Samples,Model,setting] = getData(inputFile,setting)

path = ['Inputs/',inputFile,'/model.txt'];
Model = load(path);

Model = [Model,zeros(size(Model,1),2)];

VV = (Model(2:end,2)-Model(1:end-1,2))./(Model(2:end,1)-Model(1:end-1,1));
Model(2:end-1,4) = (VV(2:end)+VV(1:end-1))/2;
Model(1,4) = Model(2,4);
Model(end,4) = Model(end-1,4);


VV = (Model(2:end,3)-Model(1:end-1,3))./(Model(2:end,1)-Model(1:end-1,1));
Model(2:end-1,5) = (VV(2:end)+VV(1:end-1))/2;
Model(1,5) = Model(2,5);
Model(end,5) = Model(end-1,5);

BEF = [Model(1,1)-0.01,-10000,0.001,0,0];
AFT = [Model(end,1)+0.01,10000,0.001,0,0];
Model = [BEF;Model;AFT];

BEF = [-10,-10000,0.001,0,0];
AFT = [40,10000,0.001,0,0];
Model = [BEF;Model;AFT];


path = ['Inputs/',inputFile,'/Data/*.txt'];
list = dir(path);

P = length(list);


data = struct('name',cell(P,1),'T',cell(P,1),'Y',cell(P,1),'R',cell(P,1));

MIN = inf;
MAX = -inf;
for p = 1:P
    data(p).name = list(p).name(1:end-4);
    
    path = ['Inputs/',inputFile,'/Data/',list(p).name];
    AA = load(path);
    data(p).T = AA(:,1);
    data(p).Y = AA(:,2);
    
    N = size(data(p).T,1);
    
    data(p).R = ones(N,1);
    
    
    MIN = min(min(data(p).T),MIN);
    MAX = max(max(data(p).T),MAX);
end

if isnan(setting.mean_T)
    setting.mean_T = (MIN+MAX)/2;
end

if isnan(setting.stdv_T)
    setting.stdv_T = (MAX-MIN)/4;
end

if isnan(setting.st)
    setting.st = min(data.T);
end

if isnan(setting.ed)
    setting.ed = max(data.T);
end

path = 'Defaults/MEAN.txt';
AA = load(path);

for p = 1:P
    data(p).T0 = linspace(setting.st,setting.ed,setting.NT)';
    
    data(p).MEAN0 = interp1(AA(:,1),AA(:,2),data(p).T0);
    data(p).MEAN = interp1(AA(:,1),AA(:,2),data(p).T);
    
    data(p).T0 = (data(p).T0-setting.mean_T)/setting.stdv_T;
    data(p).T = (data(p).T-setting.mean_T)/setting.stdv_T;
end

for p = 1:P
    data(p).SCALE = 1;
    data(p).BIAS = mean((data(p).Y-0.044)/0.033);
    data(p).BIAS = data(p).BIAS - mean(data(p).MEAN);
end


param = struct('GAMMA',cell(1,1));

param.GAMMA = zeros(2,1);
param.GAMMA(1) = 2;
param.GAMMA(2) = log(99);
param.LAMBDA = 1;


Samples = struct('name',cell(P,1));
for p = 1:P
    Samples(p).name = data(p).name;
end


end