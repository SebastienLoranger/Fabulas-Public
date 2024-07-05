classdef SM_Fwrd_Raman_2sect < MBE.System.SM_Fwrd_Raman
    %SM_FWRD_RAMAN Summary of this class goes here
    %   Single-mode Forward propagation of Raman pulse in 2 sections
    %   Only one Q-wave considered in the forward direction
    %   Only Beta and alpha can be change between section. It is assumed 
    %   that overlap and kappas are constant.
    
    properties
        startSect2(1,1) = 1;    %z position of start of section 2
        fiber2 MBE.Fiber.fiber
    end

    properties (SetAccess = protected)
        %Arrays are from  Stokes to anti-Stokes
        beta2
        alpha2
        ovlpP2
        ovlpS2
        timepos;
    end
    
    methods
        

        function obj = GenerateData(obj,fiber, medium)
            obj = GenerateData@MBE.System.SM_Fwrd_Raman(obj,fiber, medium);
            wav = 2*pi* obj.c ./ obj.omg;

            obj.beta2 = obj.fiber2.GenrateBeta(medium,wav);
            obj.alpha2 = obj.fiber2.GenrateAlpha(wav);

            [s] = obj.fiber2.GenrateOverlap(wav);
            obj.ovlpP2=zeros(1,obj.Nlines,obj.NQ);
            obj.ovlpS2=zeros(1,obj.Nlines,obj.NQ);
            for k=1:obj.NQ
                %Get overlaps from Q with most common overlap list
                %(when multiple lines are assigned to a single Q)
                stest=s(obj.Qwav{k},:);
                n=zeros(size(stest,1),1);
                for i=1:size(stest,1)
                   for j=1:size(stest,1)
                      n(i)=n(i)+double(all(stest(i,:) == stest(j,:)));
                   end
                end
                [~, m]=max(n);
                obj.ovlpS2(1,:,k)=[stest(m,:) 0];
                obj.ovlpP2(1,:,k)=[0 stest(m,:)];
            end

            
            
            obj.timepos = obj.startSect2*obj.ng/obj.c;
            
            
        end

        function obj = DefineZones(obj, mesh)
            obj = DefineZones@MBE.System.SM_Fwrd_Raman(obj,mesh);

            [~,pos] = min(abs(obj.timepos-mesh.t));
            obj.beta(pos:end,:) = repmat(obj.beta2,size(obj.beta,1)-pos+1,1);
            obj.alpha(pos:end,:) = repmat(obj.alpha2,size(obj.alpha,1)-pos+1,1);
            obj.ovlpP(pos:end,:,:) = repmat(obj.ovlpP2,size(obj.ovlpP,1)-pos+1,1,1);
            obj.ovlpS(pos:end,:,:) = repmat(obj.ovlpS2,size(obj.ovlpS,1)-pos+1,1,1);
        end



        


    end
end
