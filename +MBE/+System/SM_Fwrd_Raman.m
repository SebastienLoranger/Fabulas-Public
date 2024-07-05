classdef SM_Fwrd_Raman < MBE.System.system
    %SM_FWRD_RAMAN Summary of this class goes here
    %   Raman propagation simulation with M Q-waves.
    %   System is general for any number of wave, but SINGLE MODE
    
    properties
        MaxStokes(1,1) = 5;          %Max nb of stokes orders
        selectQ cell = {};         %Each non-empty valid cell array is a Q-wave (give line number of Stokes, 1 is lowest Stokes). All extra lines will be assigned to a general Q-wave

    end

    properties (SetAccess = protected)
        %Arrays are from  Stokes to anti-Stokes
        maxDbeta;     %Maximum coherence wave dBeta (to remove excess dB to accelerate calculation)
        kappa1;        %coupling constant: Q gain, 
        kappa2;         %coupling constant: field gain
        ovlpP;         %Overlap values applied to pump depletion (or Anti-Stokes generation)
        ovlpS;          %Overlap value applied to Stokes generation (or AS pump depletion)
        gamma;
        Qwav;     %Cell of Qs containing Stokes E-field wave index 
        phaseON;      %Row: for each Q. 1: consider PM condition. 0: auto-phase-matched, ignore exponential. (backward waves have this to 0 (auto-phase-matched))
        myStokes; %Vector (1 element per wave) indicating the wave's Stokes. If no stokes, then the element refers to current wave.
        myPump; %Vector (1 element per wave) indicating the wave's Pump. If no pump, then the element refers to current wave.
        uniform = true; %Tag indicating a single uniform section
        SingleGain = false; %if true, Q-wave will be applied ONLY to selected wave (in selectQ)
        MainQ = 0; %Index of generic Q-wave. 0 if none. Used only for identification purpose.
    end
    
    methods
        

        function obj = GenerateData(obj,fiber, medium)
            

            Order=-obj.MaxStokes:obj.MaxStokes;
            obj.Nlines=length(Order);


            
            omg_p = 2*pi * obj.c / obj.wp;
            
            neff=obj.wp * fiber.GenrateBeta(medium,obj.wp) / (2*pi);

            obj.omg = omg_p .* ones(1,length(Order)) + Order .* 2*pi .* medium.fshift(obj.wp/neff);
            obj.PumpIdx=(length(Order)-1)/2+1;

            wav = 2*pi* obj.c ./ obj.omg;

            obj.beta = fiber.GenrateBeta(medium,wav);
            obj.Nwav = length(obj.beta);
            obj.dirFwr=true(1,obj.Nwav);
            obj.ng = obj.beta(obj.PumpIdx) .* obj.wp / (2*pi);

            obj.alpha = fiber.GenrateAlpha(wav);
            kappa=medium.GenrateKappa(wav);
            obj.kappa1 = [kappa(1,:) 0];
            obj.kappa2 = [kappa(2,:) 0];

            obj.gamma = 1./medium.Dephasing(obj.wp);
            obj.maxDbeta=obj.beta(obj.PumpIdx)-obj.beta(obj.PumpIdx-1);



            %define Qs
            %Assign lines to custom Q from selectQ (validation)
            a=1:(obj.Nwav-1); obj.Qwav={[]};
            for i=1:length(obj.selectQ)
                for j=1:length(obj.selectQ{i})
                    if obj.selectQ{i}(j) > 0 && obj.selectQ{i}(j) < obj.Nwav && a(obj.selectQ{i}(j)) > 0
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
                    obj.MainQ=i+1;
                end
            end
            obj.NQ=length(obj.Qwav);
            
            %SPECIFIC FOR single-mode
            %Set overlaps for each Qs
            s = fiber.GenrateOverlap(wav); 
            obj.ovlpP=zeros(1,obj.Nwav,obj.NQ);
            obj.ovlpS=zeros(1,obj.Nwav,obj.NQ);
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
                if obj.SingleGain
                    obj.ovlpS(1,obj.Qwav{k},k)=stest(m,obj.Qwav{k});
                    obj.ovlpP(1,obj.Qwav{k}+1,k)=stest(m,obj.Qwav{k});
                else
                    obj.ovlpS(1,:,k)=[stest(m,:) 0];
                    obj.ovlpP(1,:,k)=[0 stest(m,:)];
                end
            end

            obj.maxDbeta=repmat(obj.maxDbeta,obj.NQ,1);
            obj.phaseON=ones(1,obj.Nwav,obj.NQ);
            obj.myStokes=repmat([1 1:(obj.Nwav-1)],obj.NQ,1);
            obj.myPump=repmat([2:obj.Nwav obj.Nwav],obj.NQ,1);
            
            
        end

        function obj = DefineZones(obj, mesh)
            obj = DefineZones@MBE.System.system(obj,mesh);

            obj.ovlpP=repmat(obj.ovlpP(1,:,:),length(mesh.t),1,1);
            obj.ovlpS=repmat(obj.ovlpS(1,:,:),length(mesh.t),1,1);
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
                
                %ovlp sets 0 if it is not pump or not stokes.
                f = f - 1i .* obj.ovlpP(:,:,k) .* obj.kappa2(:,obj.myStokes(k,:)) .* exp(1i.*obj.phaseON(:,:,k).* (obj.beta(:,obj.myStokes(k,:))-obj.beta+obj.maxDbeta(k)).*z) .* (obj.omg./obj.omg(:,obj.myStokes(k,:))) .* Q(:,k) .* E(:,obj.myStokes(k,:))...                 %Pump depletion (or AS gen)   
                      - 1i .* obj.ovlpS(:,:,k) .* obj.kappa2(:,:) .* exp(1i.* obj.phaseON(:,:,k).* (obj.beta(:,obj.myPump(k,:))-obj.beta-obj.maxDbeta(k)).*z)                                    .* conj(Q(:,k)) .* E(:,obj.myPump(k,:));                                          %Stokes gen (or AS pump depletion)
            end

            %NOTE: it is MUCH faster to change ALL points (10x faster!!)
            %then to change only the active points in the matrix (memory
            %access issue). HAS BEEN TESTED
        end

        function Q = DeriveQ(obj,E,z,dt)
            Q=zeros(size(E,1),obj.NQ);
            for k=1:obj.NQ
                Stk=obj.Qwav{k};
                pmp=obj.myPump(k,obj.Qwav{k});
                %RK4 in time domain for coherent wave
                for it=obj.ActiveZone(1):obj.ActiveZone(2)-1
                    %Sum of interaction between each pair of Stokes-pump
%                     if obj.uniform; iit = 1;else; iit=it;end
                    kQ = sum( 1i*0.25*obj.kappa1(it,Stk).*(E(it,pmp).*conj(E(it,Stk)).*exp(1i.* obj.phaseON(it,Stk,k) .* (obj.beta(it,pmp)-obj.beta(it,Stk)-obj.maxDbeta(k)).*z)));
                    
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

        function text = textQIDsingle(obj,i)
            if i == obj.MainQ
                text = 'Q main';
            else
                text = ['Q ' num2str(obj.Qwav{i})];
            end
        end


    end
end
