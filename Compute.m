classdef Compute
    properties
        displacement
        displacement_standardized
        v
        v_temp
        d
        I
    end
    
    methods (Static)
        function obj = principal_component(scan)
            obj.displacement = reshape(scan.displacement, numel(scan.xAxis),[]);
            obj.displacement_standardized = obj.displacement - mean(obj.displacement);
            [obj.v,obj.d] = eig(cov(obj.displacement_standardized));
            [obj.d,obj.I] = sort(diag(obj.d),1,'descend');
            obj.v = obj.v(:, obj.I);
        end
    end
end

