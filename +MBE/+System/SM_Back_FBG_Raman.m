classdef SM_Back_FBG_Raman < MBE.System.SM_Back_R_Raman
    %SM_FWRD_RAMAN Summary of this class goes here
    %   Single-mode Forward AND backward propagation of Raman pulse
    %   Only Stokes for backward signal 
    %   Applies a phase-shifted grating on first back Stokes in backward
    %   Stokes list
    
    properties
        kac(1,:) = 1;     %List  of kac. must be same size as StokesList (or ActiveBack list). If scalar, applies to first item
        Zshift(1,:)=0.5;      %List  of phase shift position (fraction of active zone). must be same size as StokesList (or ActiveBack list). If scalar, applies to first item
        FBGon(1,:)=true;     %List of logical. must be same size as StokesList (or ActiveBack list). If scalar, applies to first item
        ShiftPhi(1,:)=pi;      %List  of phase shift amplitude. must be same size as StokesList (or ActiveBack list). If scalar, applies to first item
        FBGrange(1,:)=[];      %Start-end of FBG (fraction of active zone) Default: [0 1]
    end
    
    properties (Hidden = true, SetAccess = private)
        Zkac;       %Z dependence
        FBGpos
    end

    methods
        

        function obj = DefineZones(obj, mesh)
            obj = DefineZones@MBE.System.SM_Back_R_Raman(obj,mesh);
            
            %Ensure same size of vectors. By default, only 1 FBG.
            obj.FBGon = logical(obj.FBGon);
            if size(obj.FBGon,2) < length(obj.ActiveBack)
                obj.FBGon = obj.FBGon(1) & true(1,length(obj.ActiveBack)); 
                obj.FBGon(2:end) = false;
            end
            
            %By default, only repeat parameters to other FBGs.
            if size(obj.kac,2) < sum(obj.FBGon); obj.kac = obj.kac(1) .* ones(1,sum(obj.FBGon)); end
            if size(obj.Zshift,2) < sum(obj.FBGon); obj.Zshift = obj.Zshift(1) .* ones(1,sum(obj.FBGon)); end
            if size(obj.ShiftPhi,2) < sum(obj.FBGon); obj.ShiftPhi = obj.ShiftPhi(1) .* ones(1,sum(obj.FBGon)); end
            if isempty(obj.FBGrange); obj.FBGrange=[0 1]; end
            
            
            
            
            obj.FBGpos=round(obj.FBGrange.* (obj.ActiveZone(2)-obj.ActiveZone(1))) + obj.ActiveZone(1);
            posShift = round(obj.Zshift .* (obj.FBGpos(2)-obj.FBGpos(1))) + obj.FBGpos(1);
            
            obj.Zkac = (1:length(mesh.t))';
            obj.Zkac = (obj.Zkac >= obj.FBGpos(1)) .* (obj.Zkac <= obj.FBGpos(2)) .* obj.kac;
            
            for i=1:length(posShift)
                obj.Zkac(posShift(i):end,i) = obj.Zkac(posShift(i):end,i) .* exp(1i*obj.ShiftPhi(i));
            end
        end

        function f = EqSyst(obj,E,Q,z)
            %In this system, uses all forward field + extra backard field
            %Uses single forward Q-wave + extra back Q-waves

            %Solve forward waves and backward waves seperately
            f = EqSyst@MBE.System.SM_Back_R_Raman(obj,E,Q,z);
            
            Lb = 1:length(obj.ActiveBack);

            f(:,obj.ActiveBack(obj.FBGon)) = f(:,obj.ActiveBack(obj.FBGon)) - 1i.*conj(obj.Zkac) .* E(:,obj.Nlines + Lb(obj.FBGon));
            f(:,obj.Nlines + Lb(obj.FBGon)) = f(:,obj.Nlines + Lb(obj.FBGon)) - 1i.*obj.Zkac .* E(:,obj.ActiveBack(obj.FBGon));
            %pump depletion

        end




    end
end

