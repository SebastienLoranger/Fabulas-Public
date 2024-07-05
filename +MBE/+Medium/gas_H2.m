classdef gas_H2 < MBE.Medium.gas_H2Xe
    %Fiber model
    %   Simplified waveguide optimized for lowest Stokes. 
    %   Beta not defined by fiber, but by target values
    
    properties
        Pressure(1,1)=1; %Bar
    end

    methods
  
        function obj = gas_H2()
            obj.PressureXe = 0;
        end
        
        function obj = set.Pressure(obj,val)
            obj.PressureH2=val;
        end

    end
end

