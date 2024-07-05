classdef SM_Brillouin < MBE.System.SM_Fwrd_Raman
    %SM_FWRD_RAMAN Summary of this class goes here
    %   %NOT TESTED - UNDER DEVELOPMENT
    %   Brillouin propagation simulation with M Q-waves.
    %   System is general for any number of modes
    %   Default parameter definition is for SM context
    %   DEBUG maxDbeta, check gain of SBS and PhaseON (set to 0 for
    %   forward-backward gain single Stokes shift, but 1 for reusing Q at higher orders) 
    
    
    methods
        

        function obj = GenerateData(obj,fiber, medium)
            obj.SingleGain=true;
            Nwav=2*obj.MaxStokes+1;
            obj.selectQ = {1:2:Nwav,2:2:Nwav};
            obj = GenerateData@MBE.System.SM_Fwrd_Raman(obj,fiber, medium);

            obj.dirFwr = logical(mod((1:obj.Nwav),2));
            obj.beta(~obj.dirFwr) = -obj.beta(~obj.dirFwr);
            obj.maxDbeta(1)=obj.beta(obj.PumpIdx-1)-obj.beta(obj.PumpIdx);
            obj.maxDbeta(2)=-obj.maxDbeta(1);

        end






    end
end
