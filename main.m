addpath('Codes/');

CURVE = 'BAYSPLINE'; % Define the name of calibration curve.
% CURVE = 'MULLER';
% CURVE = 'HGPR';

inputFile = 'UK37_Tierney'; % Define the input folder.
% inputFile = 'UK37_Muller';
% inputFile = 'UK37_HGPR';

kernel = 'SE[SE]'; % Define the Gaussian process kernel - SE, M25, M15, OU, SE[SE].

outputFile = CI2022(inputFile,kernel);
saveFigure(outputFile,CURVE);
