classdef General_Fwrd_Raman < MBE.System.system
    %SM_FWRD_RAMAN Summary of this class goes here
    %   NOT TESTED, UNDER DEVELOPMENT
    %   Raman propagation simulation with M Q-waves.
    %   System is general for any number of modes
    %   Default parameter definition is for SM context
    
    properties
        MaxStokes(1,1) = 5;          %Max nb of stokes orders
        selectQ cell = {};         %Each non-empty valid cell array is a Q-wave (give line number of Stokes, 1 is lowest Stokes). All extra lines will be assigned to a general Q-wave

    end

    properties (SetAccess = protected)
        %Arrays are from  Stokes to anti-Stokes
        %In this general case:
        %1st dim (rows) is time
        %2nd dim (col.) is wave
        %3rd dim (page) is Q-wave
        %4th dim is P-S duet for this Q-wave.
        maxDbeta;     %Maximum coherence wave dBeta (to remove excess dB to accelerate calculation)
        Nmodes;
        wavID;
        kappa1;        %coupling constant: Q gain, 
        kappa2;         %coupling constant: field gain
        ovlpP;         %Overlap values applied to pump depletion (or Anti-Stokes generation)
        ovlpS;          %Overlap value applied to Stokes generation (or AS pump depletion)
        gamma;
        Qwav;     %Cell of Qs containing Stokes E-field wave index (1st line) and corresponding pump (2nd line), Mode of Stokes (3rd line), Mode of pump (4th line), base-line of Stokes(5th line) and base-line of Pump (6th line)
        phaseON;      %Row: for each Q. 1: consider PM condition. 0: auto-phase-matched, ignore exponential.
        myStokes;       %1 Vector (1 column per wave) per Q per mode (4 dim total) indicating the wave's Stokes. If no stokes, then the element refers to current wave.
        myPump;         %1 Vector (1 column  per wave) per Q per mode (4 dim total) indicating the wave's Pump. If no stokes, then the element refers to current wave.
        uniform = true; %Tag indicating a single uniform section

    end
    
    methods
        

        function obj = GenerateData(obj,fiber, medium)
            

            order=-obj.MaxStokes:obj.MaxStokes;
            obj.Nlines=length(order);


            
            omg_p = 2*pi * obj.c / obj.wp;
            
            obj.omg = omg_p .* ones(1,obj.Nlines) + order .* 2*pi * medium.RamanShift;
            obj.PumpIdx=(obj.Nlines-1)/2+1;

            wav = 2*pi* obj.c ./ obj.omg;

            Beta = fiber.GenrateBeta(medium,wav);
            
            obj.Nmodes=size(Beta,1);

            obj.beta=reshape(Beta.',1,numel(Beta));
            obj.wavID=[repmat(order,1,obj.Nmodes); ceil((1:(obj.Nmodes*obj.Nlines))/obj.Nmodes)];
            obj.omg = repmat(obj.omg,1,obj.Nmodes);

            obj.Nwav = length(wav).*obj.Nmodes;
            obj.dirFwr=true(1,obj.Nwav);
            obj.ng = obj.beta(obj.PumpIdx) .* obj.wp / (2*pi);

            Alpha = fiber.GenrateAlpha(wav);
            obj.alpha=reshape(Alpha.',1,numel(Alpha));

            kappa=medium.GenrateKappa(wav);
            obj.kappa1 = repmat([kappa(1,:) 0],1,obj.Nmodes);
            obj.kappa2 = repmat([kappa(2,:) 0],1,obj.Nmodes);

            obj.gamma = 1./medium.Dephasing(wav);
            obj.maxDbeta=max(obj.beta(2:end)-obj.beta(1:end-1));


            %define Qs
            %CUSTOM ASSIGN
            %Assign lines to custom Q from selectQ (validation)
            a=1:(obj.Nlines-1); obj.Qwav={[]};
            for i=1:length(obj.selectQ)
                for j=1:length(obj.selectQ{i})
                    if obj.selectQ{i}(j) > 0 && obj.selectQ{i}(j) < obj.Nlines && a(obj.selectQ{i}(j)) > 0
                        a(obj.selectQ{i}(j))=0;
                        try
                            obj.Qwav{i}(end+1) = obj.selectQ{i}(j);
                        catch
                            obj.Qwav{i} = obj.selectQ{i}(j);
                        end
                    end
                end
            end
            %Assign remaining waves to a comon Q
            if isempty(i);i=0;end
            for k=1:length(a)
                if a(k)>0
                    try
                        obj.Qwav{i+1}(end+1) = a(k);
                    catch
                        obj.Qwav{i+1} = a(k);
                    end
                end
            end
            obj.NQ=length(obj.Qwav);

            %GENERAL ASSIGNEMENT
            %Assign related pump (line +1) for 1st mode Qs
            for i=1:obj.NQ
                obj.Qwav{i}=[obj.Qwav{i};obj.Qwav{i}+1;ones(1,size(obj.Qwav{i},2));ones(1,size(obj.Qwav{i},2));obj.Qwav{i};obj.Qwav{i}+1];
            end
            %Duplicate for Nmodes
            for j=0:obj.Nmodes-1
                for i=0:obj.Nmodes-1
                    for k=1:obj.NQ
                        if j~=0 || i~=0
                            obj.Qwav{end+1}=[obj.Qwav{k}(1,:)+j*obj.Nlines;obj.Qwav{k}(1,:)+i*obj.Nlines;j+ones(1,size(obj.Qwav{k},2));i+ones(1,size(obj.Qwav{k},2));obj.Qwav{k}(5:6,:)];
                        end
                    end
                end
            end
            obj.NQ=length(obj.Qwav);
            
            s = fiber.GenrateOverlap(wav,obj.Nmodes); 
            obj.ovlpP=zeros(1,obj.Nwav,obj.NQ,obj.Nmodes);
            obj.ovlpS=zeros(1,obj.Nwav,obj.NQ,obj.Nmodes);
            %Assign pump and Stokes
            for k=1:obj.NQ
                a=0:obj.Nmodes-1;               %mode permutation
                S=[1 1:(obj.Nlines-1)];         %1st term is not used
                P=[2:obj.Nlines obj.Nlines];    %last term is not used
                for j=0:obj.Nmodes-1
                    obj.myStokes(1,:,k,1+j)=reshape(S'+obj.Nlines*circshift(a,j),1,obj.Nmodes*obj.Nlines,1) ;
                    obj.myPump(1,:,k,1+j)=reshape(P'+obj.Nlines*circshift(a,j),1,obj.Nmodes*obj.Nlines,1) ;
                    
                    %Take 1st element (col.) of Qwav{k}
                    lineS=[];lineP=[];
                    for i=1:obj.Nmodes
                        lineS= [lineS 0 s(obj.Qwav{k}(5,1),:,fiber.getModeIdx(obj.Nmodes,obj.Qwav{k}(3,1),obj.Qwav{k}(4,1)),i+j*obj.Nmodes)];
                        lineP= [lineP s(obj.Qwav{k}(5,1),:,fiber.getModeIdx(obj.Nmodes,obj.Qwav{k}(3,1),obj.Qwav{k}(4,1)),i+j*obj.Nmodes) 0];
                    end
                    obj.ovlpS(1,:,k,j+1)=lineS;
                    obj.ovlpP(1,:,k,j+1)=lineP;
                end
            end
           
            obj.phaseON=ones(1,obj.Nwav,obj.NQ);
        end

        function obj = DefineZones(obj, mesh)
            obj = DefineZones@MBE.System.system(obj,mesh);

            obj.ovlpP=repmat(obj.ovlpP(1,:,:,:),length(mesh.t),1,1,1);
            obj.ovlpS=repmat(obj.ovlpS(1,:,:,:),length(mesh.t),1,1,1);
            obj.kappa1=repmat(obj.kappa1(1,:),length(mesh.t),1);
            obj.kappa2=repmat(obj.kappa2(1,:),length(mesh.t),1);
            obj.phaseON=repmat(obj.phaseON(1,:,:),length(mesh.t),1,1);
        end


        function f = EqSyst(obj,E,Q,z)
            %Differential equation system
            %   t: (vector) time (non-derivative, independant variable)
            %   z: (scalar) space (derivative, independant variable)
            %   y: (vector) fields to derive
            %   f: (vector) derivative of fields
  
            

            %Loss (same for all)
            f =- 0.5*obj.Zactive.*obj.alpha.*E;
            %Act along each Q
            for k=1:obj.NQ
                for p=1:obj.Nmodes
                    %ovlp sets 0 if it is not pump or not stokes.
                    f = f - 1i .* obj.ovlpP(:,:,k,p) .* obj.kappa2(:,obj.myStokes(1,:,k,p)) .* exp(1i.*obj.phaseON(:,:,k).* (obj.beta(:,obj.myStokes(1,:,k,p))-obj.beta+obj.maxDbeta).*z) .* (obj.omg./obj.omg(:,obj.myStokes(1,:,k,p))) .* Q(:,k) .* E(:,obj.myStokes(1,:,k,p))...                 %Pump depletion (or AS gen)   
                          - 1i .* obj.ovlpS(:,:,k,p) .* obj.kappa2(:,:)                     .* exp(1i.* obj.phaseON(:,:,k).* (obj.beta(:,obj.myPump(1,:,k,p))-obj.beta-obj.maxDbeta).*z)                                           .* conj(Q(:,k)) .* E(:,obj.myPump(1,:,k,p));                                          %Stokes gen (or AS pump depletion)
                end
            end

            %NOTE: it is MUCH faster to change ALL points (10x faster!!)
            %then to change only the active points in the matrix (memory
            %access issue). HAS BEEN TESTED
        end

        function Q = DeriveQ(obj,E,z,dt)
            Q=zeros(size(E,1),obj.NQ);
            for k=1:obj.NQ
                Stk=obj.Qwav{k}(1,:);
                pmp=obj.Qwav{k}(2,:);
                %RK4 in time domain for coherent wave
                for it=obj.ActiveZone(1):obj.ActiveZone(2)-1
                    %Sum of interaction between each pair of Stokes-pump
%                     if obj.uniform; iit = 1;else; iit=it;end
                    kQ = sum( 1i*0.25*obj.kappa1(it,Stk).*(E(it,pmp).*conj(E(it,Stk)).*exp(1i.* obj.phaseON(it,Stk,k) .* (obj.beta(it,pmp)-obj.beta(it,Stk)-obj.maxDbeta).*z)));
                    
                    %RK4 propagation
                    K1=-obj.gamma*Q(it,k)+kQ;
                    K2=-obj.gamma*(Q(it,k)+0.5*dt*K1)+kQ;
                    K3=-obj.gamma*(Q(it,k)+0.5*dt*K2)+kQ;
                    K4=-obj.gamma*(Q(it,k)+dt*K3)+kQ;
                
                    Q(it+1,k)=Q(it,k)+(dt/6)*(K1+2*K2+2*K3+K4);
                end


            end
        end
        
        function text = textIDsingle(obj,i)
            np=(obj.Nlines-1)/2+1;
            order = np-i;
            if order == 0
                text = 'Pump';
            elseif order>0
                text = ['Stokes ' num2str(order)];
            else
                text = ['Anti-Stokes ' num2str(-order)];
            end
        end


    end
end
