% pyversion D:\Anaconda\python.exe;
clear;clc
if count(py.sys.path, 'E:/pycharm/MyCode/Psy_tuto/centerout-valid/device/') == 0
    insert(py.sys.path, int32(0), 'E:/pycharm/MyCode/Psy_tuto/centerout-valid/device/');
end


import py.data_client.NeuracleDataClient
import py.trigger_box.TriggerNeuracle


receiver = py.data_client.NeuracleDataClient();

Control = Center_Out('par_model\TT_LDA.mat');
Control.init_figs();
pause(7);

M = 0;
for h = 1:3
    for i = 1:8
        Control.Success = 0;
        Control.Tar = Control.Tarlist(h, i);
        Control.Path = [0, 0];
        Control.Pos = [0, 0];
        Control.online_figs();
        pause(1);

        tic
        for j = 1:40
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

        % Fitt's ITR
        M = M + 1
        Control.timecost = toc;
        obj.itr = log2((0.4+0.11)/0.11) / obj.timecost;
        obj.fitts_ITR = (M-1)/M * obj.fitts_ITR + obj.itr / M;

        Control.Behave(h, i).timecost = Control.timecost;
        Control.Behave(h, i).Path = Control.Path;
        Control.Behave(h, i).Success = Control.Success;
    end
end

receiver.close();


%%
currentDate = datestr(now, 'yyyymmdd');
baseFileName = ['Log_save\Center_out_' currentDate '_'];

counter = 1;
fileName = [baseFileName num2str(counter) '.mat'];

while isfile(fileName)
    counter = counter + 1;
    fileName = [baseFileName num2str(counter) '.mat'];
end

save(fileName, "Control");

