classdef gas_SBS < MBE.Medium.medium
    %Fiber model
    %   Simplified waveguide optimized for lowest Stokes. 
    %   Beta not defined by fiber, but by target values
    %   NOT TESTED. Values and equations must be verified (Gain and dephasing time)
    
    properties
        Pressure(1,1)=1;                   %Bar
        Temperature(1,1)=298;      %Temperature (K)
        Gas = 'CO2';            %1:CO2
    end

    properties (SetAccess = protected)
        A=5e-4;                   %m^3kg^-1 (related to susceptibility)
        mass=0.04401;           %kg/mol
        sheerVisc =0;       %Sheer viscosity
        bulkVisc = 14.90e-3;       %Bulk viscosity Pa s
        thermalC = 0.01662;       %Thermal conductivity W m^-1 K^-1
        Cp=846;               %Specific heat capacity (cst pressure) [J/(K kg)]
        AdiabIndex = 1;     %Adiabatic index (heat capacity ratio)
    end

   
    methods
        function obj = gas_SBS(Gas)
            obj=obj@MBE.Medium.medium();
            obj = obj.selectGas(Gas);
            obj.freqShift = 10.8e9; %Hydrogen default shift
        end

        function obj = selectGas(obj,Gas)
            if nargin == 2
                obj.Gas=Gas;
            end
            if isequal(obj.Gas, 'CO2') || isequal(obj.Gas, 'co2')
                obj.A=5e-4;%m^3kg^-1 (Yang et al.) electrostrictive dependance of the density
                obj.mass=0.04401;%kg/mol
                obj.bulkVisc = 14.90e-3;
                obj.thermalC = 0.01662;
                obj.Cp=846;
                obj.AdiabIndex = 1.289;
            elseif isequal(obj.Gas, 'air') || isequal(obj.Gas, 'AIR')
                obj.A=5e-4;%m^3kg^-1 (copied from CO2)
                obj.bulkVisc = 14.90e-3;%WRONG (copied from CO2)
                obj.thermalC = 0.01662;%WRONG (copied from CO2)
                obj.mass=0.02897;%kg/mol
                obj.Cp=1005;
                obj.AdiabIndex = 1.4;
            elseif isequal(obj.Gas, 'ar') || isequal(obj.Gas, 'Ar')
                obj.A=5e-4;%m^3kg^-1 (copied from CO2)
                obj.bulkVisc = 14.90e-3;%WRONG (copied from CO2)
                obj.thermalC = 0.01662;%WRONG (copied from CO2)
                obj.mass=0.039948;%kg/mol
                obj.Cp=520.3;
                obj.AdiabIndex = 1.667;
            elseif isequal(obj.Gas, 'H2') || isequal(obj.Gas, 'h2')
                obj.A=5e-4;%m^3kg^-1 (copied from CO2)
                obj.bulkVisc = 14.90e-3;%WRONG (copied from CO2)
                obj.thermalC = 0.01662;%WRONG (copied from CO2)
                obj.mass=0.002016;%kg/mol
                obj.Cp=14307;
                obj.AdiabIndex = 1.405;
            elseif isequal(obj.Gas, 'He') || isequal(obj.Gas, 'he')
                obj.A=5e-4;%m^3kg^-1 (copied from CO2)
                obj.bulkVisc = 14.90e-3;%WRONG (copied from CO2)
                obj.thermalC = 0.01662;%WRONG (copied from CO2)
                obj.mass=0.004003;%kg/mol
                obj.Cp=5192.6;
                obj.AdiabIndex = 1.667;
            end

            %See https://www.ohio.edu/mechanical/thermo/property_tables/gas/idealGas.html

        end

        function va = speed(obj)
            %adiabatic compressibility 
            betas=(1/(obj.Pressure*100000))/obj.AdiabIndex; %NOT EXACT, taken for ideal gas
            rho=obj.mass*obj.density()/6.02214179e23;
            va = 1/(betas*rho).^0.5;
        end

        function f = fshift(obj,wlmat)
            %Frequency shift. Wlmat is wavelength in material
            %wlmat=wl/neff
            f=2*obj.speed()./wlmat;
        end
        
        function kappa = GenrateKappa(obj, wl, ~)
            %Vector of coupling constant between different lines.
            %Each element is the coupling between current line and its pump (omg + fr)
            %Kappa is referenced on the generated Stokes
            %1st line: gain of a coherence wave
            %2nd line: gain of a E-field 
            %Input:
            %   wl: PUMP wavelength [m]
            %   backGain (false, optional): back-gain [logical] 
            

            %From Fan Yan et all 2020 (https://doi.org/10.1038/s41566-020-0676-z)
            kappa = zeros(2,length(wl)-1);
            wl=wl(1:end-1);
            rho=obj.mass*obj.density()/6.02214179e23; %Density
            gamma=obj.A*rho; %Electrostrictive constant
            omg=2*pi*obj.c./wl;
            Tau=obj.Dephasing(wl);
            
            %Agrawal (9.4)
            
            kappa(1,:) = omg .* gamma ./ (2 .* obj.ngas(wl) .* obj.c .*rho);
            kappa(2,:) = omg .* gamma ./ (2 .* obj.c.^2 .* obj.speed());

        end

        function n = ngas(obj,lambda)
            n=1+obj.Pressure.*1e-6*(21.113+12723.2./(111-(lambda*1000000).^(-2))); %%% extracted from JOSA 67, 1551 (1977).
        end



        function T2 = Dephasing(obj,wl)
            %From Fan Yan et all 2020 (https://doi.org/10.1038/s41566-020-0676-z)
            rho=obj.mass*obj.density()/6.02214179e23;
            q=obj.fshift(wl)./obj.speed(); %wavenumber of 
            Gamma=(q.^2./rho) .* ((4/3) * obj.sheerVisc + obj.bulkVisc + obj.thermalC * (obj.AdiabIndex-1) / obj.Cp);
            T2=1./Gamma;
        end

        function D = density(obj)
            D = (obj.Pressure/1.01325)*(273.15/obj.Temperature)*2.6868e25;
        end

    end
end

