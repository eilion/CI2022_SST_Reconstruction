function [inputFile] = saveResults(inputFile,data,param,Samples,setting)

path = ['Outputs/',inputFile,'_',setting.kernel,'_',setting.var];
if exist(path,'dir') == 7
    n = 0;
    DET = 1;
    while DET == 1
        n = n + 1;
        path = ['Outputs/',inputFile,'_',setting.kernel,'_',setting.var,'(',num2str(n),')'];
        if exist(path,'dir') == 0
            DET = 0;
        end
    end
end

mkdir(path);

fileID = [path,'/results.mat'];
save(fileID,'data','param','Samples','setting');

inputFile = path(9:end);


end