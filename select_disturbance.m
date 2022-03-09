function selected_inflow = select_disturbance()

clear all;
clc;


addpath 'C:\aaa._MY_FILES\aUSN\4th semester\Vi_Forecast_data\2020April_May';
dirName = 'C:\aaa._MY_FILES\aUSN\4th semester\Vi_Forecast_data\2020April_May';

files = dir(fullfile(dirName,'*.mat'));
files = {files.name};

  
%->....50 different predicted ensemble
no_ensemble = 50;
no_days = 46;
mean_each_ensemble = zeros(no_ensemble,1);
square_error = zeros(no_ensemble,1); 
%->....Array for storing error for each 50 ensemble
selected_inflow = zeros(no_days,1);

for k = 1:no_days
    
    load(files{k});
    %->....Use all the data except row number 51 as it is control signal
    data = Vi(1:end-1,:);
    
    %->....Loop for finding mean of each 50 ensemble for 13 days
    for j = 1:no_ensemble
        val = data(j,:);
        mean_each_ensemble(j,:) = sum(val)/numel(val);
    end
    
    mean_of_mean = sum(mean_each_ensemble)/numel(mean_each_ensemble);
    
    %->....Loop for comparing mean of each ensemble to mean of mean of 
    %ensemble and finding error 
    for m = 1:no_ensemble
        square_error(m) = sqrt((mean_of_mean - mean_each_ensemble(m))^2/...
                                                             no_ensemble);   
    end
    [min_val,index] = min(square_error);
    
    selected_inflow(k) = data(index,1); 
    
end
end