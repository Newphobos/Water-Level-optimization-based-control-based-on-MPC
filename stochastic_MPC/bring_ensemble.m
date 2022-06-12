function  rep_ensemble = bring_ensemble(day_no,no_ensemble,N_days)

    addpath 'C:\aaa._MY_FILES\aUSN\4th semester\Vi_Forecast_data\2020April_May_NormalData';
    dirName = 'C:\aaa._MY_FILES\aUSN\4th semester\Vi_Forecast_data\2020April_May_NormalData';

    files = dir(fullfile(dirName,'*.mat'));
    files = {files.name};

    load(files{day_no});
            %->....Use all the data except row number 51 as it is control signal
    ensem1 = Vi(1:no_ensemble,1:N_days);
       
    rep_ensemble = repelem(ensem1,1,24);
    
end