classdef fiber
    %Fiber model
    %   Simplified waveguide optimized for lowest Stokes. 
    %   Fundamental mode only.
    %   Beta not defined by fiber, but by target values
    %   Overlap between last Stokes and all other comb lines can be
    %   adjusted. It is supposed that all other overlaps are 1.
    
    properties

        Dcore = 60e-6;       %Core diameter

    end

    properties (Hidden = true)
        c=299792458;
    end

    methods(Static)
        function index = getModeIdx(Nmode,S,P)
            a=circshift(0:Nmode-1,P-1);index = a(S)*Nmode+S;
        end
    end
    
    methods
        

        function B = GenrateBeta(~,medium, wl)
            B=2*pi*medium.ngas(wl)./wl;
        end

        function alpha = GenrateAlpha(~, wl)
            alpha=zeros(1,length(wl));
        end

        function [s] = GenrateOverlap(~,wl,Nmode)
            %Overlap between current generated Stokes and its pump
            %s is 2 or 4 dimension matrix
            %dim 1 (lines): pump order of Q (overlaps with m-1 Stokes field)
            %dim 2 (col.): pump order of E-field (overlaps with k-l Stokes field)
            %
            %Note: Matrix will be: (Nlines-1 x  Nlines-1)
            %Ex: single-mode, 5 lines (example from Loranger et al.)
            %   P:      S1     P       AS1     AS2
            %   Q       
            %   S1      s2     s1      s1      s1
            %   P       s1      1       1       1
            %   AS1     s1      1       1       1
            %   AS2     s1      1       1       1
            %
            %For multimode:
            %dim 3: Mode combination for Q
            %dim 4: Mode combination for E-field
            %
            %For dim 3 and 4:
            %Index position of element Stokes s with pump p:
            %a=circshift(0:Nmode-1,p-1);index = a(s)*Nmode+s;
            %Ex for 3 modes:
            %Mode Index: 1       2       3       4       5       6       7       8       9
            %s-p        1-1     2-2     3-3     1-3     2-1     3-2     1-2     2-3     3-1
            %function getModeIdx(Nmode,S,P) can be used
            %Note: in each mode index position, there are Nlines
            if nargin<3
                Nmode=1;
            end
            Nstokes=length(wl);
            s=ones((Nstokes-1),(Nstokes-1),Nmode.^2,Nmode.^2);
            

        end

        function A = Area(obj)
            A = pi.*(obj.Dcore/2).^2;
        end

        
    end
end

