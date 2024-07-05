classdef Fields
    % Properties of fiber.
    % Type .help for details
    
    properties

        %Field parameters
        Pump_energy(1,1) = 100e-6;   %J
        Pulse_width(1,1) = 3e-9;     %s
        noisefloor(1,:) = 5;         %V/m
        Impedance(1,1) = 377;        %Only for display W conversion. Not critical to calculation
        
        ErrTol(1,1) = 1e-3;          %Field error tolerance (V/m)
        textID
        textQID
    end
    
    properties (SetAccess = private)       
        E;       %complete Field (time, lines, position)
        Q;
        z;
        noiseMat;
    end

    properties (Hidden = true)
        c=299792458
        eps0 = 8.85418781e-12
        h = 6.62607015e-34
        iz = 1      %Current counter for z save
        area
        Qz;      %Current coherence wave at current position z
        Ez;      %Current field at current position z
    end
    
    methods(Static)
        function help()
            fprintf('\nDescription:\n')
        end
    end
    
    methods
        

        function obj = InitFields(obj,mesh, PumpIdx, Nwav, NQ)
            %Initialize field
            %Temporal field E is in [V/m]
            %Spatial field F is unitless []
            %Power defined as P = (1/2eta) |E|^2 Aeff = 
            if length(obj.noisefloor) < Nwav
                obj.noisefloor = obj.noisefloor(1) .* ones(1,Nwav);
            end
            
            pump_width_I=obj.Pulse_width/sqrt(2*log(2));  % temporal width of the Gaussian pump [s]
            P_peak=obj.Pump_energy/(sqrt(0.5*pi)*pump_width_I); % pump peak power [W] 
            
            
            obj.E=zeros(length(mesh.t),Nwav, length(mesh.z));
            obj.Q=zeros(length(mesh.t),NQ, length(mesh.z));
            
            obj.Qz=zeros(length(mesh.t),NQ);
            obj.Ez = obj.noisefloor.*ones(length(mesh.t),Nwav);
%             obj.Ez(:,1:Nlines-NQ+1) = obj.noisefloor(1:Nlines-NQ+1).*exp(-((mesh.t-mesh.t_start)./pump_width_I).^2);

            obj.Ez(:,PumpIdx) = obj.Ez(:,PumpIdx)+sqrt(2*P_peak/(obj.eps0*obj.c*obj.area)).*exp(-((mesh.t-mesh.t_start)./pump_width_I).^2);% pump pulse electric field envelope [V/m]

            obj.E(:,:,1) = obj.Ez;
            obj.iz=2;
            obj.z=mesh.z;
            obj.noiseMat = repmat(obj.noisefloor,size(obj.Ez,1),1);

        end
        
        function obj = MinNoise(obj)
           select = abs(obj.Ez) < obj.noiseMat;
           obj.Ez(select) = obj.noiseMat(select);
        end

        function I = Intensity(obj)
            %Result in W
            %each row (dim 1): position in time
            %each col (dim 2): position in z
            %each page (dim 3): field
            I = 0.5*permute(abs(obj.E).^2,[1 3 2])/obj.Impedance;
        end

        function Qr = ReshapeQ(obj)
            %Result in W
            %each row (dim 1): position in time
            %each col (dim 2): position in z
            %each page (dim 3): field
            Qr = 0.5*permute(abs(obj.Q),[1 3 2]);
        end
        
        function phi = Phase(obj)
            %Result in W
            %each row (dim 1): position in time
            %each col (dim 2): position in z
            %each page (dim 3): field
            phi = permute(angle(obj.E),[1 3 2]);
        end
        

        function J = Energy(obj,dt,reverse)
            %dt: time step. Result in J 
            %reverse: integrate on slices instead of time
            %each row (dim 1): position in Z
            %each col (dim 2): field
            if nargin==2
                reverse = false;
            end
            if reverse
                J = 0.5*sum(abs(obj.E).^2,3)*dt/obj.Impedance;
            else
                J = 0.5*sum(permute(abs(obj.E).^2,[3 2 1]),3)*dt/obj.Impedance;
            end
        end

        function N = Photon(obj,dt,freq,reverse)
            %dt: stime step
            %freq: vector (row) of requencies
            %result in number of photons.
            %each row (dim 1): position in Z
            %each col (dim 2): field
            if nargin==3
                reverse = false;
            end
            N = obj.Energy(dt,reverse)./(obj.h * freq);
        end
        
        function obj = KeepE(obj,zcur)
            if obj.iz<= size(obj.E,3)
                
                %Interpolate to correct for out-of-grid point
                obj.E(:,:,obj.iz) = obj.Ez - (zcur-obj.z(obj.iz)) .* (obj.Ez-obj.E(:,:,obj.iz-1))./(zcur-obj.z(obj.iz-1));
                obj.Q(:,:,obj.iz) = obj.Qz - (zcur-obj.z(obj.iz)) .* (obj.Qz-obj.Q(:,:,obj.iz-1))./(zcur-obj.z(obj.iz-1));
%                 obj.E(:,:,obj.iz) = obj.Ez;
%                 obj.Q(:,obj.iz) = obj.Qz(:,1);

                obj.iz = obj.iz+1;
            end
        end
        
        function obj = update(obj,z)
            obj.E=cat(3,obj.E,zeros(size(obj.E,1),size(obj.E,2),length(z)-size(obj.E,3)));
            obj.Q=cat(3,obj.Q,zeros(size(obj.Q,1),size(obj.Q,2),length(z)-size(obj.Q,3)));
            obj.z=z;
        end

    end
    
    
end