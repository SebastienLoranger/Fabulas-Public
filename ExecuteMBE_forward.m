%% Initiate model
clear(); close all
%Load system 
model=MBE.Model;
model.saveFile='';
model.syst = MBE.System.SM_Fwrd_Raman;      %Select system
model.fiber = MBE.Fiber.Ideal_S;            %Select waveguide
model.medium = MBE.Medium.gas_H2Xe;           %Select medium
model.syst.MaxStokes = 2;                   %Max nb of stokes orders
model.syst.selectQ={[1]};                   %Each non-empty valid cell array is an extra Q-wave (give line number of Stokes, 1 is lowest Stokes). All extra lines will be assigned to a general Q-wave

%% Enter parametesr
%Pump
GenFreq=3;  %Thz
model.syst.wp = 299792458/(GenFreq*10^12+2*124.57e12);  %Pump wavelength
model.field.Pump_energy = 500e-6;                       %J
model.field.Pulse_width = 3e-9;                         %s
model.field.noisefloor = [30 30 30 30 30];                             %V/m
model.field.Impedance = 0.5;                          %Only for display W conversion. Not critical to calculation


%Waveguide
model.fiber.loss = [0 0 0 0 0];               %Loss (dB/m). Scalar: uniform to all modes; Vector: specific loss for each mode
model.fiber.dBetaS = 75;             % Phase-matching offset (m^-1) from ideal case
model.fiber.dBetaAS = [1000 2000];         % Phase mismatched for AS. Scalar: uniform to all AS. Vector: specific to each AS line.
model.fiber.overlap = [0.35;0.5767];            % Overlap integrals [s2; s1] 
model.fiber.thick = 500e-9;         %Core-wall thickness
model.fiber.Dcore = 80e-6;          %Core diameter

% model.syst.fiber2=model.fiber;
% model.syst.fiber2.dBetaS = 1000; 
% model.syst.startSect2 = 0.4; %m

%Gas
model.medium.PressureH2=70;            %Bar
model.medium.PressureXe=0;            %Bar
model.medium.Temperature=298;       %Temperature (K)

%simulation parameters
model.adaptive = true;                      %Adapts dz to acceptable error
model.mesh.RelativeWindow = false;          %True: centered on pulse. Faster but less accurate (entrance dynamics not taken into account). Cannot be used for counter propagating pulses or if presence of interface.
model.mesh.Nt = 2^12;                       % Number of time points of pulse (for non-relative window, total number of points is increased depending on length)
model.mesh.Width_t = 12e-9;                 % Window width of pulse (s) (for non-relative window, total time is increased depending on length)
        
model.mesh.L_fib = 3;                           %Spatial window (m) (for non-relative only): this will significantly increase time window (and points)
model.mesh.L_sim = 3;                       %Simulation lenght (m) (stop simulation after this propagation)
        
model.mesh.dz_start = 1e-6;                 %Starting step size. If Adaptive if FALSE, the this is definitive step size.
model.mesh.dz_floor=10e-9;                    % minimum dz step in adaptive calculation (m). If reaches minimum, then calculation is less accurate (error too high).
model.mesh.Nsave = 300;                     %Number of z-slices to save
model.field.ErrTol = 1e-3;                  %Field error tolerance (V/m) for adaptive calculation



model.PrepareSim(); %Initialises simulation
model.run(11,[1 2 3]); %Run simulation. Note: can be called again to continue if Z increased or after paused. Arguments are fore display purpose only.
model.displayEz(100,{'b','r','k'})
model.displayEt(0,[1 1 1], 1)
model.displayQ(10, 1)

%Note: for old versions, do: Model.field.textID = Model.syst.textID();
        

%Note: to continue propagation, do:
% model.mesh.L_sim = L_new; %Where L_new is new propagation length (longer then initial L_sim)
% model.run();