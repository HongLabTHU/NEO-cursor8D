% pyversion D:\Anaconda\python.exe;
clear;clc
if count(py.sys.path, 'E:/pycharm/MyCode/Psy_tuto/centerout-valid/device/') == 0
    insert(py.sys.path, int32(0), 'E:/pycharm/MyCode/Psy_tuto/centerout-valid/device/');
end


import py.data_client.NeuracleDataClient
import py.trigger_box.TriggerNeuracle


receiver = py.data_client.NeuracleDataClient();
pause(3);


Control = Web_grid('par_model\TT_LDA.mat');


for i = 1:64
    Control.Success = 0;
    Control.Tar = Control.Tarlist(i);
    Control.Path = Control.Pos;
    Control.Tar_old = Control.Pos;
    Control.online_figs();

    tic
    for j = 1:70
        a = toc;
        data = cell(receiver.get_trial_data());
        X = double(data{3}.flatten().tolist());
        X = reshape(X, 1000, 8)' / 1e-6;
        X = X - repmat(mean(X), size(X, 1), 1);
        [spectra, fb] = Online_psd(X, 1000, 100, 256);

        Control = Control.online_task(spectra);

        if Control.Success
            break;
        end

        pause(0.2 - (toc - a) - 0.004);
        [j , toc - a]
    end

    Control.Behave(i).timecost = toc;
    Control.Behave(i).Path = Control.Path;
    Control.Behave(i).Success = Control.Success;
    Control.hit_figs(Control.Behave(i).timecost);
end


receiver.close();


%%
currentDate = datestr(now, 'yyyymmdd');
baseFileName = ['Log_save\Web_grid_' currentDate '_'];

counter = 1;
fileName = [baseFileName num2str(counter) '.mat'];

while isfile(fileName)
    counter = counter + 1;
    fileName = [baseFileName num2str(counter) '.mat'];
end

save(fileName, "Control");


