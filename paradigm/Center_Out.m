classdef Center_Out
    properties
        Tarlist        
        Tar
        Pos
        classer

        Behave
        timecost
        Success
        Path

        AX
    end
    
    methods
        function obj = Center_Out(classerpath)
            obj.classer = load(classerpath);
            obj.Tarlist = [randperm(8); randperm(8); randperm(8)];
            obj.Tar = 1;
            obj.Pos = [0, 0];

            obj.Success = 0;
            obj.timecost = NaN;
            obj.Path = [0, 0];

            obj.Behave(3, 8).timecost = obj.timecost;
            obj.Behave(3, 8).Path = obj.Path;
            obj.Behave(3, 8).Success = obj.Success;

            close all;
            figure("Position", [2769,-233,1524,920], "Color", [1 1 1])
            obj.AX = axes('Position', [0 0 1 1]);
        end
        
        function obj = online_task(obj, spectra)

            tarpos = 0.4*[cos((obj.Tar-1)*pi/4) sin((obj.Tar-1)*pi/4)];
            V2 = obj.velocity_model(spectra);
            obj.Pos = obj.Pos + 0.025 * V2;

            if norm(obj.Pos) > 0.41 
                obj.Pos = obj.Pos / norm(obj.Pos) * 0.41;
            end

            obj.Path = [obj.Path; obj.Pos];
            if norm(obj.Pos - tarpos) < 0.11
                obj.Success = 1;
            end
            obj.online_figs();
        end

        function online_figs(obj)
            cla(obj.AX);
            r = 0.06;
            for i = setdiff(1:8, obj.Tar)
                tarpos = 0.4*[cos((i-1)*pi/4) sin((i-1)*pi/4)];
                pos = [tarpos(1)-r, tarpos(2)-r, 2*r, 2*r];
                rectangle(obj.AX, 'Position', pos, 'Curvature', [1 1], 'LineWidth', 1);
            end

            tarpos = 0.4*[cos((obj.Tar-1)*pi/4) sin((obj.Tar-1)*pi/4)];
            pos = [tarpos(1)-r, tarpos(2)-r, 2*r, 2*r];
            rectangle(obj.AX, 'Position', pos, 'Curvature', [1 1], 'FaceColor', [0 0.4470 0.7410], 'LineWidth', 1);
 
            r = 0.05;
            pos = [obj.Pos(1)-r, obj.Pos(2)-r, 2*r, 2*r];
            rectangle(obj.AX, 'Position', pos, 'Curvature', [1 1], 'FaceColor', [0.8500 0.3250 0.0980], 'LineWidth', 1);

            axis(obj.AX, [-0.48, 0.48, -0.48, 0.48]);
            axis(obj.AX, "equal");
            axis(obj.AX, "off");
            hold(obj.AX, 'off');
        end

        function init_figs(obj)
            cla(obj.AX);
            text(0.2, 0.05, '光标控制准备开始', 'FontSize', 50);
            axis(obj.AX, [-0.48, 0.48, -0.48, 0.48]);
            axis(obj.AX, "equal");
            axis(obj.AX, "off");
        end

        function V2 = velocity_model(obj, spectra)
            X = log10(spectra);
            X = (X(:)' - obj.classer.meanX) ./ obj.classer.stdX;
            [pred, ~, ~, posterior] = predict(obj.classer.classer, X);

            V = posterior(:, 2:end) * obj.classer.V_turning;
            V1 = (1-posterior(:, 1).^0.25) .* V;
            V2 = [real(V1) imag(V1)];
        end
    end
end