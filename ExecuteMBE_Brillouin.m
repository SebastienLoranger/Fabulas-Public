clear(); %close all
%NOT TESTED - UNDER DEVELOPMENT
%Load system 
model=MBE.Model;
model.saveFile='';
model.syst = MBE.System.SM_Brillouin;      %Select system
model.fiber = MBE.Fiber.Ideal_S;            %Select waveguide
model.medium = MBE.Medium.gas_SBS('CO2');           %Select medium
model.syst.MaxStokes = 2;                   %Max nb of stokes orders

%Pump

model.syst.wp = 1.55e-6;  %Pump wavelength
model.field.Pump_energy = 500e-6;                       %J
model.field.Pulse_width = 3e-9;                         %s
model.field.noisefloor = [50 50 50 50 50];                             %V/m
model.field.Impedance = 0.5;                          %Only for display W conversion. Not critical to calculation


%Waveguide
model.fiber.loss = [0 0 0 0 0];               %Loss (dB/m). Scalar: uniform to all modes; Vector: specific loss for each mode
model.fiber.dBetaS = 0;             % Phase-matching offset (m^-1) from ideal case
model.fiber.dBetaAS = [0 0];         % Phase mismatched for AS. Scalar: uniform to all AS. Vector: specific to each AS line.
model.fiber.overlap = 1;            % Non-linear overlap (acoustic wave and pump)
% model.fiber.overlap = [1; 1];              %Very ideal case 
model.fiber.thick = 500e-9;         %Core-wall thickness
model.fiber.Dcore = 80e-6;          %Core diameter

% model.syst.fiber2=model.fiber;
% model.syst.fiber2.dBetaS = 1000; 
% model.syst.startSect2 = 0.4; %m

%Gas
model.medium.Pressure=70;            %Bar
model.medium.Temperature=298;       %Temperature (K)

%simulation parameters
model.adaptive = true;
model.mesh.RelativeWindow = false;
model.mesh.Nt = 2^12;                       % Number of time points (of pulse only)
model.mesh.Width_t = 12e-9;                 % Window width (s) (of pulse only)
        
model.mesh.L_fib = 10;                           %Simulation length (m)
model.mesh.L_sim = 12;
% model.mesh.Nt = 2^12;                       % Number of time points
% model.mesh.Width_t = 20e-9;                 % Window width (s)
        
% model.mesh.L_fib = 1;                           %Simulation length (m)
model.mesh.dz_start = 1e-3;                 %Starting step size
model.mesh.dz_floor=10e-9;                    % minimum dz step tolerable (m)
model.mesh.Nsave = 300;  

model.field.ErrTol = 1e-3;                  %Field error tolerance (V/m)


model.PrepareSim();
model.run(11,[1 2 3]);
model.displayEz(100,{'b','r','k'})
model.displayEt(0,[1 1 1], 1)
model.displayQ(10, 1)


        

%Note: to continue propagation, do:
% model.mesh.L_sim = L_new; %Where L_new is new propagation length (longer then initial L_sim)
% model.run();