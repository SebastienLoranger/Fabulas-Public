classdef Model < handle

    properties
        syst MBE.System.system = MBE.System.SM_Fwrd_Raman;
        medium MBE.Medium.medium = MBE.Medium.gas_H2;
        fiber MBE.Fiber.fiber = MBE.Fiber.Ideal_S;
        mesh MBE.Data.Mesh = MBE.Data.Mesh;
        field MBE.Data.Fields = MBE.Data.Fields;
        adaptive logical = true;
        saveFile char
    end
    properties (Hidden = true, SetAccess = private)
        c=299792458;
        eps0 = 8.85418781e-12;
        iit;             %Time iteration in moving window
        relSpeed;       %Relative speed factor for back-wave. 1 if non relative, 2 if relative
    end

    properties (Hidden = true)
        cmap{};           %color map for display    
    end

    properties (SetAccess = private)
        Bad;            %Triggers to True if dz reaches acceptable floor
        z;              %Current propagation position

    end
   
    methods(Static)
       %% documentation
        function help
            fprintf('\nProperties of model:\n')

        end
    end
   
    methods
        %Constructor
        function obj = Model()
            data = load('+MBE\+Data\cmap_biegert.mat');
            obj.cmap{1}=data.cmap_biegert;
            obj.cmap{4}=[linspace(1,0.2,100)' linspace(1,0,100)' linspace(1,1,100)'];
            obj.cmap{2}=[linspace(1,0,100)' linspace(1,0.5,100)' linspace(1,0,100)'];
            obj.cmap{3}=[linspace(1,1,100)' linspace(1,0.2,100)' linspace(1,0,100)'];
        end
       
       %Prepare grid
        function PrepareSim(obj)

            if isempty(obj.mesh.RelativeWindow)
                obj.mesh.RelativeWindow = obj.adaptive;
            end
            if isempty(obj.mesh.ForceIntStep)
                obj.mesh.ForceIntStep = ~obj.adaptive;
            end

            obj.syst = obj.syst.GenerateData(obj.fiber, obj.medium);

            obj.mesh = obj.mesh.run(obj.syst.ng);
            obj.syst = obj.syst.DefineZones(obj.mesh);
            

            obj.field.area = obj.fiber.Area();
            obj.field = obj.field.InitFields(obj.mesh, obj.syst.PumpIdx, obj.syst.Nwav,obj.syst.NQ);
            obj.Bad=false;
            obj.z=obj.mesh.z(1);
            obj.iit = 1;
            obj.field.textID = obj.syst.textID();
            obj.field.textQID = obj.syst.textQID();
            obj.relSpeed = 1 + double(obj.mesh.RelativeWindow);
        end
        
        
        function run(obj, fig,selectE,selectQ)
            %run()
            %run(fig)
            %run(fig, selectE)
            %run(fig, selectE, selectQ)
            %fig: figure number (optional) to display progress
            %selectE: vector of fields to display (optional). Default: all
            %fields
            %selectQ: vector of Q-wave to display (optional). Default: all
            %waves
            Emax=max(max(max(abs(obj.field.E))));
            if nargin>1
                if nargin < 4
                    selectQ = 1:size(obj.field.Qz,2);
                    if nargin<3
                        selectE = 1:size(obj.field.Ez,2);
                    end
                end
                figure(fig);legendText=[obj.field.textID(selectE) obj.field.textQID(selectQ)];
                h=plot(obj.mesh.t,20*log10(abs(obj.field.Ez(:,selectE(1))./Emax)));hold on
                for i=2:length(selectE)
                    h(end+1)=plot(obj.mesh.t,20*log10(abs(obj.field.Ez(:,selectE(i))./Emax)));
                end
                for i=1:length(selectQ)
                    h(end+1)=plot(obj.mesh.t,20*log10(abs(obj.field.Qz(:,selectQ(i)))/0.009));
                end
                legend(legendText);
                hold off;
            else 
                h=[];
            end
            
            my_clock=cputime;
            if obj.z>obj.mesh.z(1)
                fprintf('\n Resuming simulation...');
                if obj.mesh.L_sim>obj.mesh.z(end)
                   obj.mesh=obj.mesh.update(); 
                   obj.field=obj.field.update(obj.mesh.z);
                end
            else
                fprintf('\n Starting simulation...');
            end
            
            
            while obj.z < obj.mesh.z(end)
                %Propagate
                if obj.adaptive
                    obj.IterateAdaptive();
                else
                    obj.IterateFixed();
                end
                
                %Displace time window for z-dependant (non-relative case) or fo any backward propagating waves
                while obj.z > (obj.mesh.t(obj.iit) - obj.mesh.t(1))*obj.relSpeed*obj.c/obj.syst.ng && obj.iit < length(obj.mesh.t)
                    obj.iit=obj.iit+1;
                    obj.field.Ez = obj.syst.MoveField(obj.field.Ez,obj.mesh.RelativeWindow);
                end
                
                %Break on divergence
                if any(any(isnan(obj.field.Ez)))
                    fprintf('\n Diverging field! Abort simulation');
                    break
                end

                %Save data slice
                if obj.z > obj.mesh.z(obj.field.iz)
                    obj.field=obj.field.KeepE(obj.z);
                    fprintf('\nz = %f m        percentage completed : %f %%        dz = %e',obj.z, obj.z*100/obj.mesh.z(end), obj.mesh.dz);
                    if obj.Bad; fprintf('WARNING: Reached min dz! field error beyond tolerance'); end
                    if ~isempty(h)
                        for i=1:length(selectE)
                            h(i).YData=20*log10(abs(obj.field.Ez(:,selectE(i))./Emax));
                        end
                        for j=1:length(selectQ)
                            h(i+j).YData=20*log10(abs(obj.field.Qz(:,selectQ(j)))./0.009);
                        end
                        pause(0.01)
                    end
                end
            end

            fprintf('\n Finished after %f min\n\n',(cputime-my_clock)/60);
            
            if ~isempty(obj.saveFile)
                obj.Save();
            end

        end
        
        function IterateFixed(obj)
            %dz propagation iteration
            %WARNING! This method supposes same group velocity (same group
            %index) for each pulse. Valid only for long pulses at
            %near-phase-match (very close neff, hence similar ng)

            %Refine RK4 propagation to compare error
            obj.field.Qz = obj.syst.DeriveQ(obj.field.Ez,obj.z, obj.mesh.dt);
            obj.field.Ez = obj.rk4(obj.field.Ez, obj.field.Qz, obj.mesh.dz);
            obj.z=obj.z+obj.mesh.dz;
            obj.field = obj.field.MinNoise();

        end
        
        function IterateAdaptive(obj)
            %dz propagation iteration

            %Refine RK4 propagation to compare error
            Q1 = obj.syst.DeriveQ(obj.field.Ez,obj.z, obj.mesh.dt);
            E1 = obj.rk4(obj.field.Ez, Q1, obj.mesh.dz);
            E2 = obj.rk4(obj.field.Ez, Q1, obj.mesh.dz/2);
            Q2 = obj.syst.DeriveQ(E2, obj.z, obj.mesh.dt);
            E3 = obj.rk4(E2, Q2, obj.mesh.dz/2);

            %Calculate error
            Err = max(sqrt(sum(abs(E3-E1).^2)./sum(abs(E3).^2)));
            
            %adaptive step (adjust dz according to error)
            if Err > 2*obj.field.ErrTol && obj.mesh.dz > obj.mesh.dz_floor
                obj.mesh.dz = obj.mesh.dz/2;
            else
                if Err < 0.5*obj.field.ErrTol && obj.mesh.dz < obj.mesh.max_dz
                    obj.mesh.dz = 1.148698354997035 * obj.mesh.dz;
                elseif Err > obj.field.ErrTol && Err <= 2*obj.field.ErrTol
                    obj.mesh.dz = obj.mesh.dz/1.148698354997035;
                elseif obj.mesh.dz <= obj.mesh.dz_floor
                    obj.Bad=true;
                end

                %Advance to next step
                obj.z=obj.z+obj.mesh.dz;
                obj.field.Ez = E3;
                obj.field.Qz = Q2;
                obj.field = obj.field.MinNoise();
            end
            
        end

        function E = rk4(obj, E, Q, dz)
           k1 = obj.syst.EqSyst(E, Q, obj.z);
           k2 = obj.syst.EqSyst(E + dz*k1/2, Q, obj.z + dz/2); 
           k3 = obj.syst.EqSyst(E + dz*k2/2, Q, obj.z + dz/2);
           k4 = obj.syst.EqSyst(E + dz*k3, Q, obj.z + dz);
           E = E + (dz/6) * (k1+2*(k2+k3) + k4);
        end


        function Save(Model,varargin)
            if nargin >= 2
                file = varargin{1};
            else
                file = Model.saveFile;
            end

            try 
                save(file, 'Model');
            catch
                fprintf('\nWarning: Unable to save, invalid file path\n')
            end
        end


        function displayEz(obj, varargin)
            %Displays normalized photon count along z-axis
            %
            %displayEz() 
            %displayEz(fig) 
            %displayEz(fig, cline)
            %displayEz(fig, cline, reverse)
            %
            %fig (scalar):  figure number. Default is new figure.
            %cline (cell-vector): each element is a color property
            %       for the corresponding wave. Empty string indicates wave 
            %       will not be displayed. Size of cell-vector can be less
            %       then number of waves, in which case the not-present
            %       waves will not be displayed.
            %       Ex: cline = {'k','','b'}
            %           k: black   --: dashed lines
            %           b: blue    : : dotted lines
            %           g: green   . : single dot points
            %           r: red     o : single 'o' points
            %reverse (bool): Reverses axis. true:energy is integrated on
            %                 slices instead of time.
            reverse = false;
            if nargin>1
                figure(varargin{1});
                if nargin>2
                    cline = varargin{2};
                    if nargin>3
                        reverse = varargin{3};
                    end
                else
                    cline ={};
                end
            else
                figure();
            end
            if reverse
                dt=(obj.mesh.z(2)-obj.mesh.z(1))/obj.c;
            else
                dt=obj.mesh.dt;
            end
            Phot = obj.field.Photon(dt,obj.syst.omg/(2*pi),reverse);
            maxP = max(sum(Phot,2));
            text={};
            if isempty(cline)
                for i=1:size(Phot,2)
                    if reverse
                        plot(obj.mesh.t*obj.c,100*Phot(:,i)./maxP,'LineWidth',2);hold on;
                    else
                        plot(obj.mesh.z,100*Phot(:,i)./maxP,'LineWidth',2);hold on;
                    end
                    text{i}=obj.field.textID{i};
                end
            else
                for i=1:min([size(Phot,2),length(cline)])
                    if ~isempty(cline{i})
                        if reverse
                            plot(obj.mesh.t*obj.c,100*Phot(:,i)./maxP, cline{i}, 'LineWidth',2);hold on;
                        else
                            plot(obj.mesh.z,100*Phot(:,i)./maxP, cline{i}, 'LineWidth',2);hold on;
                        end
                        text{end+1}=obj.field.textID{i};
                    end
                end
            end
            hold off;
            legend(text);
            xlabel('Distance [m]','Fontsize',15); 
            ylabel('Normalized photon numbers','Fontsize',20); 
            set(gca,'Fontsize',15);
        end

        function displayEt(obj, varargin)
            %Displays E-waves in time-space 2D color plot
            %displayEt() 
            %displayEt(fig)  
            %displayEt(fig, display) 
            %displayEt(fig, display, reverse) 
            %displayEt(fig, display, reverse, spaceTime) : 
            %
            %fig (scalar):  displays in figure number "fig". if empty, 
            %               displays in new figure.
            %display (vector bool): Displays selected fields. Item
            %       size above the number of fields will be ignored. First 
            %       item is first field.
            %reverse (bool): Reverses axis. X-axis is time-data
            %       (converted to space) and Y-axis is space slices 
            %       (converted to time)
            %spaceTime (bool): for relative-window only (will be ignore 
            %       if non-relative). Displays in absolute time vs space 
            %       (tilted figure) instead of relative default time-frame.
            

            spaceTime = false;reverse =false;
            if nargin>1
                fig = varargin{1};
                if nargin>2
                    display = varargin{2};
                    if nargin>3
                        reverse = varargin{3};
                        if nargin>4
                            spaceTime=varargin{4};
                        end
                    end
                else
                    display = ones(1,obj.syst.Nwav);
                end
            else
                fig = [];
            end
            
            Int = obj.field.Intensity();
            xtext='Time (ns)';
            if obj.mesh.RelativeWindow && spaceTime
                t_disp=obj.mesh.t(1):obj.mesh.t(2)-obj.mesh.t(1):(obj.mesh.z(end)/obj.c+obj.mesh.t(end));
                Int_zt = zeros(length(t_disp),size(Int,2),size(Int,3));
                for j=1:obj.mesh.Nsave
                    [~, it] = min(abs(t_disp-(obj.mesh.z(j)/obj.c-obj.mesh.t(end))));
                    Int_zt(1+it:length(obj.mesh.t)+it,j,:)=Int(:,j,:);
                end
                Int = Int_zt;
            else 
                t_disp=obj.mesh.t;
                if obj.mesh.RelativeWindow; xtext='Relative time (ns)'; end
            end

            for i=1:min([size(Int,3),length(display)])
                if display(i)
                    if isempty(fig)
                        h=figure();
                    else
                        h=figure(fig+i);
                    end
                    if reverse
                        imagesc(t_disp.*obj.c,10^9*obj.mesh.z./obj.c,Int(:,:,i).');colormap(obj.cmap{1});
                    else
                        imagesc(obj.mesh.z,t_disp,Int(:,:,i));colormap(obj.cmap{1});
                    end
                    title(obj.field.textID{i},'Fontsize',20);xlabel('Distance (m)','Fontsize',15); ylabel(xtext,'Fontsize',15); set(gca,'Fontsize',15);
                    set(h, 'DefaultFigurePosition', [50+i*100 300 405 332]); %%% set the position/size of the figure on the screen 
                end
            end

        end

        function displayQ(obj, varargin)
            %Displays Q-waves in time-space 2D color plot
            %displayQ()
            %displayQ(fig) 
            %displayQ(fig, reverse)
            %displayQ(fig, reverse, spaceTime, xlim)
            %
            %fig (scalar):  displays in figure number "fig". if empty, 
            %               displays in new figure.
            %reverse (bool): Reverses axis. X-axis is time-data
            %       (converted to space) and Y-axis is space slices 
            %       (converted to time)
            %spaceTime (bool): for relative-window only (will be ignore 
            %       if non-relative). Displays in absolute time vs space 
            %       (tilted figure) instead of relative default time-frame.
            %xlimit (vector): vector of x-limit.

            spaceTime=false;reverse =false;xlimit=[];
            if nargin>1
                fig = varargin{1};
                if nargin>2
                    reverse = varargin{2};
                    if nargin>3
                        spaceTime=varargin{3};
                        if nargin>4
                            xlimit=varargin{4};
                        end
                    end
                end
            else
                fig = [];
            end


            if isempty(fig)
                h=figure();
            else
                h=figure(fig);
            end

            Q = obj.field.ReshapeQ();
            xtext='Time (ns)';
            if obj.mesh.RelativeWindow && spaceTime
                t_disp=obj.mesh.t(1):obj.mesh.t(2)-obj.mesh.t(1):(obj.mesh.z(end)/obj.c+obj.mesh.t(end));
                Q_zt = zeros(length(t_disp),size(Q,2),size(Q,3));
                for j=1:obj.mesh.Nsave
                    [~, it] = min(abs(t_disp-(obj.mesh.z(j)/obj.c-obj.mesh.t(end))));
                    Q_zt(1+it:length(obj.mesh.t)+it,j,:)=Q(:,j,:);
                end
                Q = Q_zt;
            else 
                if obj.mesh.RelativeWindow; xtext='Relative time (ns)'; end
                t_disp=obj.mesh.t;
            end
            set(h, 'Position', [400 200 500 400]); %%% set the position/size of the figure on the screen 
            maxQ=max(max(max(Q)));vector=[];
            if size(Q,3)<4

                for k=1:size(Q,3)
                    ax{k} = axes;
                    if reverse
                        him{k}=imagesc(ax{k},t_disp.*obj.c,10^9*obj.mesh.z./obj.c,10*log10(Q(:,:,size(Q,3)-k+1)./maxQ).'); hold on
                    else
                        him{k}=imagesc(ax{k},obj.mesh.z,t_disp,10*log10(Q(:,:,size(Q,3)-k+1))); hold on
                    end
                   title('Coherence wave','Fontsize',20);xlabel('Distance (m)','Fontsize',15); ylabel(xtext,'Fontsize',15); set(gca,'Fontsize',15);
                    ax{k}.CLim=[-20 0];
                    T=(10*log10(Q(:,:,size(Q,3)-k+1)./maxQ)>-20) .* ((10*log10(Q(:,:,size(Q,3)-k+1)./maxQ)+20)/20);
                    T(isnan(T))=0;
                    set(him{k},'AlphaData',T.');
                    vector=[vector ax{k}];
%                     axis square; 
                     ax{k}.Position=[0.15,0.18,0.65,0.70];
                     if ~isempty(xlimit)
                        xlim(ax{k},xlimit)
                     end
                end
                linkaxes(vector) ;
                ax{k}.Visible='off';
                ax{k}.XTick = []; 
                ax{k}.YTick = [];
                for k=1:size(Q,3)
                    colormap(ax{k},obj.cmap{1+k});
                    cbar=colorbar(ax{k},'location','east');
                    cbar.Position=[0.82+(k-1)*0.07    0.18    0.0427    0.7];
                    cbar.Ticks=[-20 -15 -10 -5 0];
                    cbar.TickLabels=[];
                end
%                 set(vector,'Position',[.17 .11 .685 .815]);
                hold off

            else
                if reverse
                    imagesc(t_disp.*obj.c,10^9*obj.mesh.z./obj.c,sum(Q,3)); colormap(obj.cmap{1});
                else
                    imagesc(obj.mesh.z,t_disp,sum(Q,3));colormap(obj.cmap{1});
                end
                title('Coherence wave','Fontsize',20);xlabel('Distance (m)','Fontsize',15); ylabel(xtext,'Fontsize',15); set(gca,'Fontsize',15);
                
            end
            
            

        end

        function Display(obj)
            obj.displayEz(100,{'b','r','k'},1)
            obj.displayEt(0,[1 1 1], 1)
            obj.displayQ(10, 1)
        end
    

    end

end