classdef gas_H2Xe < MBE.Medium.medium
    %Fiber model
    %   Simplified waveguide optimized for lowest Stokes. 
    %   Beta not defined by fiber, but by target values
    
    properties
        PressureH2(1,1)=1;                   %Bar
        PressureXe(1,1)=1;                   %Bar
        
        Temperature(1,1)=298;      %Temperature (K)
        backGainRatio(1,1) = 0.4;  %Backward gain ratio (empirical value based on experiment)
    end

   
    methods
        function obj = gas_H2Xe()
            obj=obj@MBE.Medium.medium();
            obj.freqShift = 124.57e12; %Hydrogen default shift
        end

        function kappa = GenrateKappa(obj, wl, backGain)
            %Vector of coupling constant between different lines.
            %Each element is the coupling between current line and its pump (omg + fr)
            %Kappa is referenced on the generated Stokes
            %1st line: gain of a coherence wave
            %2nd line: gain of a E-field 
            %Input:
            %   wl: PUMP wavelength [m]
            %   backGain (false, optional): back-gain [logical] 
            freqs=obj.c./wl;
            gain = obj.GainH2(freqs(2:end)); %gain at pump wavelength
            if nargin==3 && backGain
                gain=gain*obj.backGainRatio;
            end
            gamma=1./obj.Dephasing();

            kappa = zeros(2,length(wl)-1);
            kappa(1,:) = -sqrt((2*obj.c.^2.*obj.eps0.^2.*gain.*gamma)./(obj.density()*obj.hbar*2*pi*freqs(1:end-1))); 
            kappa(2,:) = obj.density()*obj.hbar*2*pi*freqs(1:end-1).*abs(kappa(1,:))/(2*obj.c*obj.eps0); 

        end

        function n = ngas(obj,lambda)
            n=1+obj.PressureH2.*1e-6*(21.113+12723.2./(111-(lambda*1000000).^(-2))); %%% extracted from JOSA 67, 1551 (1977).
        end

        function gain = GainH2(obj, frep) % pressure H2, pressure buffer gas, temperature, pump frequency, all in SI units
            rho_ama = (obj.PressureH2/1.01325)*(273.15/obj.Temperature); % number density in amagat
            K_B = 0.658; % Boltzmann population of J=1 level at 298 K
            cvak = 299792458; % speed of light in vacuum
            nu_MHz = 1/(pi*obj.Dephasing())*1e-6; % linewidth in MHz
            nu_p = frep/cvak/100; % pump frequency in inverse cm (wavenumber)
            gain_cm = (9.37e6).*(52*rho_ama./nu_MHz).*(K_B/0.658).*(nu_p-4155.25).*(7.19e9-nu_p.^2).^(-2); % steady-state gain in cm/W
            gain = gain_cm/100; % gain in m/W
        end

        function T2 = Dephasing(obj,~)
            if obj.PressureXe==0
                rho_ama = (obj.PressureH2/1.01325)*(273.15/obj.Temperature); % number density in amagat
                nu_MHz = (309/rho_ama)*(obj.Temperature/298)^0.92+(51.8+0.152*(obj.Temperature-298)+(4.85e-4)*(obj.Temperature-298)^2)*rho_ama; % linewidth (FWHM) in MHz
                nu = nu_MHz*1e6; % linewidth (FWHM) in Hz
                T2 = 1/(pi*nu);  % dephasing time in s
            else
                rho_ama1 = (obj.PressureH2/1.01325)*(273.15/obj.Temperature); % number density in amagat
                rho_ama2 = (obj.PressureXe/1.01325)*(273.15/obj.Temperature); % number density in amagat
                rho_1=rho_ama1; rho_2=rho_ama2;
                rho_t=rho_1+rho_2; % Total pressure as the sum of the hydrogen pressure (1) and buffer gas pressure (2)
                Dif=1./(rho_1/(1.42*rho_t)+rho_2/(0.5*rho_t)); % difusion coefficient for a binary gas mixture
                Acoef=(((2*pi*4155)^2)*Dif/pi)*1e-6; % A coefficient
                Bcoef1=51.8; Bcoef2=414; % B coefficients for H2 (1) and Xe (2)
                nu_MHz_H2_Xe=Acoef/rho_t+(Bcoef1*rho_1+Bcoef2*rho_2);
                nu = nu_MHz_H2_Xe*1e6; % linewidth (FWHM) in Hz
                T2 = 1/(pi*nu);  % dephasing time in s
            end
        end

        function D = density(obj)
            D = (obj.PressureH2/1.01325)*(273.15/obj.Temperature)*2.6868e25;
        end

    end
end

