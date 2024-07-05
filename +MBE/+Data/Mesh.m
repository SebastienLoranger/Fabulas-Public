classdef Mesh
    % Properties of fiber.
    % Type .help for details
    
    properties
        %Sim parameters
        Nt(1,1) = 2^12;                       % Number of time points of pulse region 
        Width_t(1,1) = 20e-9;                 % Window width of pulse (s)
        
        L_fib(1,1) = 1;                           %Active zone length (m)                    
        L_sim(1,:) = [];                     %Simulation length
        dz_start(1,1) = 1e-6;                 %Starting step size
        dz_floor(1,1)=1e-12;                    % minimum dz step tolerable (m)
        Nsave(1,1) = 200;                    %Number of slices (z) to save

        RelativeWindow(1,:) logical = [];           %True: centered on pulse. Faster but less accurate (entrance dynamics not taken into account). Cannot be used for counter propagating pulses or if presence of interface.

        ForceIntStep(1,:) logical = [];         %Force dz selection to integer number of dt
    end
    
    properties (Hidden = true)
        dz; % when properties defining kappas were last modified
        c=299792458;
        t_start=0;
    end

    properties (SetAccess = private)
        t;
        z;
        max_dz;
        dt;
        ActiveZone;     %vector of 2 elements (start-end) of gain medium zone (in position #)
    end
    
    methods(Static)
        function help()
            fprintf('\nDescription:\n')
        end
    end
    
    methods
        function obj = run(obj,varargin)
            %Generate grid
            %run()
            %run(ng)
            %   ng: group index (default is 1)
            %


            if nargin >1
                ng=varargin{1};
            else
                ng=1;
            end
            
            %Set time mesh
            if ~isempty(obj.RelativeWindow) && obj.RelativeWindow(1)
                obj.t = ((0:1:obj.Nt-1)' - (obj.Nt/2)) * obj.Width_t / obj.Nt;
                obj.dt = obj.t(2) - obj.t(1);
                obj.ActiveZone = [1 length(obj.t)];
                obj.t_start=0;
            else
                obj.dt = obj.Width_t / obj.Nt;
                obj.t = (-obj.Width_t : obj.dt : obj.L_fib*ng/obj.c)' ;
                [~, obj.ActiveZone(1)] = min(abs(obj.t));
                [~, obj.ActiveZone(2)] = min(abs(obj.t - obj.L_fib*ng/obj.c));
                obj.t_start=-obj.Width_t/2;
            end
            
            %Set initial dz
            if ~isempty(obj.ForceIntStep) && obj.ForceIntStep(1)
                obj.dz = ceil(obj.dz_start / (obj.dt*obj.c/ng)) .* obj.dt*obj.c/ng;
            else
                obj.dz = obj.dz_start;
            end
            
            %Set slicing elements for saving
            if isempty(obj.L_sim)
                obj.L_sim(1)=obj.L_fib;
            end
            obj.z=linspace(0,obj.L_sim,obj.Nsave);
            obj.max_dz=(obj.z(2)-obj.z(1))/1.15;
        end
        
        function obj=update(obj,varargin)
            %Update z slicing for longer simulation length
            %update()
            %   Uses current mesh.L_sim
            %update(L_new)
            %   L_new: New length. Will update mesh.L_sim.
            
            if nargin >1
                obj.L_sim=varargin{1};
            end
            
            Dz=obj.z(2)-obj.z(1);
            obj.z=[obj.z (obj.z(end)+Dz):Dz:obj.L_sim];
            obj.Nsave=length(obj.z);
        end
    end
    
    
end