function [outputFile] = CI2022(inputFile,kernel)

disp('-----------------------------------------------------------------');
disp('## This is a multivariate Gaussian process state-space model.');
disp(['#  Input file is ',inputFile,'.']);
disp(['#  Kernel is ',kernel,'.']);

% Load setting:
setting = getSetting(inputFile,kernel);

% Load data:
[data,param,Samples,Model,setting] = getData(inputFile,setting);

% Get individual samples:
disp('#  Sampling point-wise SST samples...');
[data,Samples] = getSamples_INDIV(data,Samples,Model);
disp('   Done.');

% Learn the transition model:
disp('#  Learning parameters...');
[data,param] = learnParameters(data,param,setting);
disp('   Done.');

% Draw HMC samples:
disp('#  Sampling GPST SST samples...');
[data,Samples] = getSamples(data,param,Samples,Model,setting);
[data,Samples] = getSamples_Final(data,param,Samples,setting);
disp('   Done.');

% Restore Inputs:
data = restoreVariables(data,setting);

outputFile = saveResults(inputFile,data,param,Samples,setting);
% saveFigure(outputFile);
disp(['## Results are saved in the folder Outputs/',outputFile,'/.']);
disp('-----------------------------------------------------------------');


end
