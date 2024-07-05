%% Initiate model
clear(); %close all
model=MBE.Model;
model.saveFile='FBG-DFB10m-1-500uJ-6ns-Loss.mat';
model.syst = MBE.System.SM_Back_FBG_Raman;
model.fiber = MBE.Fiber.Ideal_S;            %Select waveguide
model.medium = MBE.Medium.gas_H2Xe;           %Select medium
model.adaptive = true;
model.mesh.RelativeWindow = false;
model.syst.MaxStokes = 2;                   %Max nb of stokes orders
model.syst.StokesList=[2];                % List of backward Stokes orders (ex: [2 1] is 2nd and 1st Stokes, negative numbers are AntiStokes. cannot be highest anti-Stokes)
%% Enter parametesr
%Pump
GenFreq=3;  %Thz
model.syst.wp = 299792458/(GenFreq*10^12+2*124.57e12);  %Pump wavelength
model.field.Pump_energy = 500e-6;                       %J
model.field.Pulse_width = 6e-9;                         %s
model.field.noisefloor = 1*[10 50 50 10 10 0];                             %V/m
model.field.Impedance = 0.5;                          %Only for display W conversion. Not critical to calculation


%Waveguide
model.fiber.loss = [10 0 0 0 0 10];               %Loss (dB/m). Scalar: uniform to all modes; Vector: specific loss for each mode
model.fiber.dBetaS = 0;             % Phase-matching offset (m^-1) from ideal case
model.fiber.dBetaAS = [1000 8000];         % Phase mismatched for AS. Scalar: uniform to all AS. Vector: specific to each AS line.
model.fiber.overlap = [1 1];            %Overlap between pump and lowest Stokes. 
model.fiber.thick = 500e-9;         %Core-wall thickness
model.fiber.Dcore = 80e-6;          %Core diameter

%FBG
model.syst.kac = 10;     %List of kac. must be same size as StokesList (or ActiveBack list). If scalar, applies to first item
model.syst.Zshift=0.1122;      %List of phase shift position (fraction of active zone). must be same size as StokesList (or ActiveBack list). If scalar, applies to first item
model.syst.ShiftPhi=pi;      %List of phase shift amplitude. must be same size as StokesList (or ActiveBack list). If scalar, applies to first item
model.syst.FBGrange=[0 1];

%Gas
model.medium.PressureH2=70;             %Bar
model.medium.PressureXe=0;              %Bar
model.medium.Temperature=298;           %Temperature (K)
model.medium.backGainRatio = 0.4;         %Backward gain ratio (empirical value based on experiment)

%simulation parameters
model.mesh.Nt = 2^14;                       % Number of time points (of pulse only)
model.mesh.Width_t = 20e-9;                 % Window width (s) (of pulse only)
model.mesh.L_fib = 3;                       %Length of active region
model.mesh.L_sim = 8;                       %Simulation length (m) (propagation)
model.mesh.dz_start = 1e-6;                 %Starting step size
model.mesh.dz_floor=1e-12;                    % minimum dz step tolerable (m)
model.mesh.Nsave = 200;  
model.field.ErrTol = 1e-3;                  %Field error tolerance (V/m)

%% Run simulation
model.PrepareSim();
model.run(1,[3 2 1 6]);
model.displayEz(100,{'b','r','k','','','g--'})
model.displayEt(0,[1 1 1 0 0 1], 1)
model.displayQ(10, 1)
        