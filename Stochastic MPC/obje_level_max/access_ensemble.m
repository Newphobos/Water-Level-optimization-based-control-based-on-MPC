function all_ensemble = access_ensemble()
    



    clear all;
    clc;

    Np = 312;
    addpath 'C:\aaa._MY_FILES\aUSN\4th semester\Vi_Forecast_data\2020April_May';
    dirName = 'C:\aaa._MY_FILES\aUSN\4th semester\Vi_Forecast_data\2020April_May';

    files = dir(fullfile(dirName,'*.mat'));
    files = {files.name};

    no_ensemble = 50;
    no_days = 30;
    
    all_ensemble = zeros(50*Np,no_days);
    
    for k = 1:no_days
    
            load(files{k});
            %->....Use all the data except row number 51 as it is control signal
            data = Vi(1:end-1,:);

           
            for j = 1:no_ensemble
                ensem = data(j,:);

                ensem_repeat = repelem(ensem,24);

                all_ensemble((j-1)*Np + 1:j*Np,k) = ensem_repeat;
            end

    end
    
end

