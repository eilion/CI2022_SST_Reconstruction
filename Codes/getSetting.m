function [setting] = getSetting(inputFile,kernel)

setting = struct('max_iters',cell(1,1));

path = 'Defaults/setting.txt';
fileID = fopen(path);
INFO = textscan(fileID,'%s %s');
fclose(fileID);

setting.nSamples = str2double(INFO{2}{strcmp(INFO{1},'number_of_samples:')==1});
setting.max_iters = str2double(INFO{2}{strcmp(INFO{1},'max_iters:')==1});
setting.kernel = kernel;
setting.mean_T = str2double(INFO{2}{strcmp(INFO{1},'mean_T:')==1});
setting.stdv_T = str2double(INFO{2}{strcmp(INFO{1},'stdv_T:')==1});
setting.st = str2double(INFO{2}{strcmp(INFO{1},'query_age_start:')==1});
setting.ed = str2double(INFO{2}{strcmp(INFO{1},'query_age_end:')==1});
setting.NT = str2double(INFO{2}{strcmp(INFO{1},'number_of_query_ages:')==1});


path = ['Inputs/',inputFile,'/setting.txt'];
if exist(path,'file') == 2
    fileID = fopen(path);
    INFO = textscan(fileID,'%s %s');
    fclose(fileID);
    
    if sum(strcmp(INFO{1},'number_of_samples:')==1) == 1
        setting.nSamples = str2double(INFO{2}{strcmp(INFO{1},'number_of_samples:')==1});
    end
    
    if sum(strcmp(INFO{1},'max_iters:')==1) == 1
        setting.max_iters = str2double(INFO{2}{strcmp(INFO{1},'max_iters:')==1});
    end
    
    if sum(strcmp(INFO{1},'mean_T:')==1) == 1
        setting.mean_T = str2double(INFO{2}{strcmp(INFO{1},'mean_T:')==1});
    end
    
    if sum(strcmp(INFO{1},'stdv_T:')==1) == 1
        setting.stdv_T = str2double(INFO{2}{strcmp(INFO{1},'stdv_T:')==1});
    end
    
    if sum(strcmp(INFO{1},'query_age_start:')==1) == 1
        setting.st = str2double(INFO{2}{strcmp(INFO{1},'query_age_start:')==1});
    end
    
    if sum(strcmp(INFO{1},'query_age_end:')==1) == 1
        setting.ed = str2double(INFO{2}{strcmp(INFO{1},'query_age_end:')==1});
    end
    
    if sum(strcmp(INFO{1},'number_of_query_ages:')==1) == 1
        setting.NT = str2double(INFO{2}{strcmp(INFO{1},'number_of_query_ages:')==1});
    end
end


end