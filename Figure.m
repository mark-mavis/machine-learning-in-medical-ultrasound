classdef Figure
    methods (Static)
        function bMode(scan, figure_title)
            bModeTemp = log10(mean(scan.bMode,3));
            figure('Name', figure_title, 'NumberTitle', 'on');
            clf;
            pcolor(scan.xAxis,scan.zAxis,bModeTemp);
            axis ij;
            axis equal;
            axis tight;
            shading flat;
            caxis(prctile(bModeTemp(:),[1,99]));
            colormap(gray);
        end

        function scan_masks(scan, figure_title)
            % ctLabelled: 3d matrix made up of all masks
            ctLabelled = scan.skullMaskThick;
            ctLabelled(:,:,2) = scan.bloodMaskThick;
            ctLabelled(:,:,3) = scan.ventMaskThick;
            
            figure('Name', figure_title, 'NumberTitle', 'on');
            clf;
            H = pcolor(scan.xAxis,scan.zAxis, ctLabelled(:,:,1));
            axis ij;
            axis equal;
            axis tight;
            shading flat;

            % Assigning ctLabelled 3d matrix made up of CT scan probabilities
            % to the H. CData
            H.CData = ctLabelled;
        end

        function principle_components_figure(scan, figure_title, pc_obj)
            figure('Name', strcat('principle components', ' | ', figure_title), 'NumberTitle', 'on');
            clf;
            
            for loop1 = 1:9
                subplot(3,3,loop1);
                plot(scan.timeAxis,pc_obj.v(:,loop1));
                title(sprintf('PC %d (%0.1f%% variance)',loop1, pc_obj.d(loop1)/sum(pc_obj.d)*100));
            end
        end

        function pcScore(scan, figure_title, pc_obj)
            figure('Name', strcat('pcScore', ' | ', figure_title), 'NumberTitle', 'on');
            clf;
            for loop1 = 1:9
                subplot(3,3,loop1);
                pcScore = pc_obj.displacement_standardized*pc_obj.v(:,loop1);
                pcScore = reshape(pcScore,size(scan.xAxis));
                pcolor(scan.xAxis,scan.zAxis,pcScore);
                axis ij;axis equal;axis tight;shading flat;
                caxis(prctile(pcScore(:),[1,99]));
                colormap(parula);
            end
        end
    end
end

