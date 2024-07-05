classdef medium
    %Base medium class
    
    properties
        freqShift = 124.57e12;       %Raman frequency shift (Hz)
    end

    properties (Hidden = true)
        c=299792458;
        eps0 = 8.85418781e-12;
        hbar=1.054571726e-34;
    end
    
    methods
        function f = fshift(obj,~)
            %Frequency shift. Wlmat is wavelength in material
            %wlmat=wl/neff
            f=obj.freqShift;
        end

        function kappa = GenrateKappa(~, wl, ~)
            kappa = zeros(2,length(wl)-1);
        end
        
        function n = ngas(~,~)
            n=1; %%% extracted from JOSA 67, 1551 (1977).
        end

        function T2 = Dephasing(~,~)
           T2 = 1;  % dephasing time in s
        end
    end
end

