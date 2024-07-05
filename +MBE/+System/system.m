classdef system 
    %SM_FWRD_RAMAN Summary of this class goes here
    %   Single-mode Forward propagation of Raman pulse
    %
    
    properties
        wp(1,1) = 1030e-9;           %Pump wavelength    
    end
    
    properties (Hidden = true, SetAccess = protected)
        c=299792458;
        ActiveZone=[];
        Zactive;       %Z dependence
    end

    properties (SetAccess = protected)
        %Arrays are from  Stokes to anti-Stokes
        omg;          % angular frequencies
        ng;
        Nwav;       %Number of waves in the system (all directions)
        Nlines;     %Number of lines in the system (stokes orders + pump)
        NQ;         %Number of coherent waves
        alpha;        %loss (vector, size = nb of line)
        beta;   %propagation constants (vector, size = nb of line). 
        dirFwr;     %Vector of logical (1 col. for each eave). 1: forward. 0: Backward 
        PumpIdx; %Index number of pump wave
    end
    
    methods
        

        function obj = GenerateData(obj,fiber, medium) %Dummy
            obj.Nwav=1;
            obj.NQ=1;
            omg_p = 2*pi * obj.c / obj.wp;
            obj.omg = omg_p; wav = 2*pi* obj.c ./ obj.omg;
            obj.beta = fiber.GenrateBeta(medium,wav);
            obj.ng = obj.beta .* obj.wp / (2*pi);
            obj.Nlines=length(obj.beta);
            obj.PumpIdx=1; 
        end

        function obj = DefineZones(obj, mesh)
            obj.ActiveZone = mesh.ActiveZone;
            obj.Zactive = (1:length(mesh.t))';
            obj.Zactive = (obj.Zactive >= obj.ActiveZone(1)) .* (obj.Zactive <= obj.ActiveZone(2));
            obj.beta=repmat(obj.beta(1,:),length(mesh.t),1);
            obj.alpha=repmat(obj.alpha(1,:),length(mesh.t),1);
        end


        function f = EqSyst(~,E,~,~) %Dummy
            %Differential equation system
            %   t: (vector) time (non-derivative, independant variable)
            %   z: (scalar) space (derivative, independant variable)
            %   y: (vector) fields to derive
            %   f: (vector) derivative of fields
            f = zeros(size(E,1),size(E,2));
            
        end

        
        function Q = DeriveQ(~,E,~,~) %Dummy
            Q=ones(size(E,1),1);
        end
        
        function E = MoveField(obj,E,relative)
            if ~relative
                E(:,obj.dirFwr) = circshift(E(:,obj.dirFwr),1,1);
                E(1,obj.dirFwr) = 0.*E(1,obj.dirFwr);
            end
            E(:,~obj.dirFwr) = circshift(E(:,~obj.dirFwr),-1,1);
            E(end,~obj.dirFwr) = 0 .* E(end,~obj.dirFwr);
        end
        
        function text = textID(obj,i)
            if nargin > 1
                text = obj.textIDsingle(i);
            else
                text=cell(1,obj.Nwav);
                for i=1:obj.Nwav
                    text{i} = obj.textIDsingle(i);
                end
            end
        end

        function text = textQID(obj,i)
            if nargin > 1
                text = obj.textQIDsingle(i);
            else
                text=cell(1,obj.NQ);
                for i=1:obj.NQ
                    text{i} = obj.textQIDsingle(i);
                end
            end
        end
        

        function text = textIDsingle(~,i) %Dummy
            text = ['Wave ' num2str(i)];
        end

        function text = textQIDsingle(~,i) %Dummy
            text = ['Q ' num2str(i)];
        end


    end
end
