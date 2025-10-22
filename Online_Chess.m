classdef Online_Chess < handle
    properties
        Tarlist        
        Tar
        Tar_old
        Pos
        classer
        Mesh1
        Mesh2

        Behave
        timecost
        Success
        Path

        AX
    end
    
    methods
        function obj = Online_Chess(classerpath)
            obj.classer = load(classerpath);
            obj.Tarlist = randperm(64);
            obj.Tar = 1;
            obj.Pos = [0, 0];
            obj.Tar_old = [0, 0];

            obj.Success = 0;
            obj.timecost = NaN;
            obj.Path = [0, 0];
            obj.Mesh1 = meshgrid(-0.35:0.1:0.35);
            obj.Mesh2 = meshgrid(-0.35:0.1:0.35)';

            obj.Behave(64).timecost = obj.timecost;
            obj.Behave(64).Path = obj.Path;
            obj.Behave(64).Success = obj.Success;

            close all;
            figure("Position", [769,-233,1524,920], "Color", [1 1 1])
            obj.AX = axes('Position', [0 0 1 1]);
        end
        
        function obj = online_task(obj, spectra)
            % 光标控制
            tarpos = [obj.Mesh1(obj.Tar) obj.Mesh2(obj.Tar)];
            V2 = obj.velocity_model(spectra);
            obj.Pos = obj.Pos + 0.025 * V2;

            obj.Pos(1) = max(min(obj.Pos(1), 0.42), -0.42);
            obj.Pos(2) = max(min(obj.Pos(2), 0.42), -0.42);

            obj.Path = [obj.Path; obj.Pos];
            obj.online_figs();
            if norm(obj.Pos - tarpos) < 0.08
                obj.Success = 1;
            end
            
        end

        function hit_figs(obj, timecost) % figures after hit
            if Control.Success
                tarpos = [obj.Mesh1(obj.Tar) obj.Mesh2(obj.Tar)];
                r = 0.045;
                pos = [tarpos(1)-r, tarpos(2)-r, 2*r, 2*r];
                rectangle(obj.AX, 'Position', pos, 'Curvature', [1 1], 'FaceColor', [0.9290 0.6940 0.1250], 'LineWidth', 1);
                D = norm([obj.Mesh1(obj.Tar) obj.Mesh2(obj.Tar)] - obj.Tar_old);
                Fitts = log2((D+0.085)/0.085) / timecost;
                text(obj.AX, 0.5, 0.4, ['Fitts ITR ' num2str(Fitts, '%.2f') 'bpm'], 'FontSize', 20);
                pause(0.3);
            end
        end

        function online_figs(obj)
            cla(obj.AX);
            r = 0.045;
            for i = 1:9
                CCC = 0.85*(abs(i-5)<4) * ones(1, 3);
                plot(obj.AX, [0.1*i-0.5 0.1*i-0.5], [-0.40 0.40], 'Color', CCC, 'LineWidth', 1);
                hold(obj.AX, 'on');
                plot(obj.AX, [-0.40 0.40], [0.1*i-0.5 0.1*i-0.5], 'Color', CCC, 'LineWidth', 1);
            end

            tarpos = [obj.Mesh1(obj.Tar) obj.Mesh2(obj.Tar)];
            pos = [tarpos(1)-r, tarpos(2)-r, 2*r, 2*r];
            rectangle(obj.AX, 'Position', pos, 'Curvature', [1 1], 'FaceColor', [0 0.4470 0.7410], 'LineWidth', 1);
 
            r = 0.04;
            pos = [obj.Pos(1)-r, obj.Pos(2)-r, 2*r, 2*r];
            rectangle(obj.AX, 'Position', pos, 'Curvature', [1 1], 'FaceColor', [0.8500 0.3250 0.0980], 'LineWidth', 1);

            axis(obj.AX, [-0.45, 0.45, -0.45, 0.45]);
            axis(obj.AX, "equal");
            axis(obj.AX, "off");
            hold(obj.AX, 'off');
        end

        function init_figs(obj)
            cla(obj.AX);
            text(0.2, 0.05, '光标控制准备开始', 'FontSize', 50);
            axis(obj.AX, [-0.45, 0.45, -0.45, 0.45]);
            axis(obj.AX, "equal");
            axis(obj.AX, "off");
        end

        function V2 = velocity_model(obj, spectra)
            % 速度控制
            X = log10(spectra);
            X = (X(:)' - obj.classer.meanX) ./ obj.classer.stdX;

            [pred, ~, ~, posterior] = predict(obj.classer.classer, X);
            % [pred, ~, posterior, ~, ~] = classify(X, obj.classer.X, obj.classer.ylab, 'linear');

            V_turning = [-sqrt(2) + 1j*sqrt(2);
                -sqrt(2) - 1j*sqrt(2);
                sqrt(2)  + 1j*sqrt(2);
                sqrt(2)  - 1j*sqrt(2);
                -2       + 1j*0      ;
                -0       + 1j*2      ;
                -0       - 1j*2      ;
                2        + 1j*0      ];


            V = posterior(:, 2:end) * V_turning;
            V2 = (1-posterior(:, 1).^0.25) .* V;
            V2 = [real(V1) imag(V1)];
        end
    end
end