classdef SM_Back_R_Raman < MBE.System.SM_Fwrd_Raman
    %SM_FWRD_RAMAN Summary of this class goes here
    %   Single-mode Forward AND backward propagation of Raman pulse
    %   Only Stokes for backward signal 
    
    properties
        StokesList(1,:)=[];          % List of backward Stokes orders (ex: [2 1] is 2nd and 1st Stokes, negative numbers are AntiStokes. cannot be highest anti-Stokes)

        Reflect(1,:) = [];          %Reflectivity of end-of-fiber
    end
    


    properties (Hidden = true)
        ActiveBack;
    end

    methods
        

        function obj = GenerateData(obj,fiber, medium)
            
            

            obj = GenerateData@MBE.System.SM_Fwrd_Raman(obj,fiber, medium);
            wav = 2*pi* obj.c ./ obj.omg;

            %Backward overlap
            
            %Add Q-waves
            k=1;
            for i=1:length(obj.StokesList)
                a= (obj.Nlines-1)/2 + 1 - obj.StokesList(i);
                if a >= 1 && a < obj.Nlines
                    obj.ActiveBack(end+1) = a;
                    obj.Qwav{obj.NQ+k}=obj.Nwav+k;
                    k=k+1;
                end
            end

            obj.beta = [obj.beta obj.beta(obj.ActiveBack)];
            obj.alpha= [obj.alpha obj.alpha(obj.ActiveBack)];
            obj.dirFwr=[obj.dirFwr false(1,length(obj.ActiveBack))];
            %Update all other arrays
            Nback=length(obj.ActiveBack);
            obj.Nwav=obj.Nlines+Nback;
            obj.ovlpP=cat(3,[obj.ovlpP zeros(1,Nback,obj.NQ)],zeros(1,obj.Nwav,Nback));
            obj.ovlpS=cat(3,[obj.ovlpS zeros(1,Nback,obj.NQ)],zeros(1,obj.Nwav,Nback));

            obj.omg=[obj.omg obj.omg(obj.ActiveBack)];

            kappaB=medium.GenrateKappa(wav,true);
            obj.kappa1 = [obj.kappa1 kappaB(1,obj.ActiveBack)];
            obj.kappa2 = [obj.kappa2 kappaB(2,obj.ActiveBack)];
            
            overlapB=diag(fiber.GenrateOverlap(wav)).';

            obj.myStokes=[obj.myStokes ones(obj.NQ,Nback);ones(Nback,obj.Nlines+Nback)];
            obj.myPump=[obj.myPump ones(obj.NQ,Nback);ones(Nback,obj.Nlines+Nback)];
            for i=1:Nback
                obj.ovlpS(1,obj.Nlines+i,obj.NQ+i)=overlapB(obj.ActiveBack(i));
                obj.ovlpP(1,obj.ActiveBack(i)+1,obj.NQ+i)=overlapB(obj.ActiveBack(i));

                obj.myStokes(obj.NQ+i,obj.ActiveBack(i)+1) = obj.Nlines+i;
                obj.myPump(obj.NQ+i,obj.Nlines+i) = obj.ActiveBack(i)+1;
            end

            obj.phaseON=cat(3,[obj.phaseON zeros(1,Nback,obj.NQ)], zeros(1,obj.Nlines+Nback,Nback));
            obj.NQ=length(obj.Qwav);
            
            obj.maxDbeta=[obj.maxDbeta;zeros(Nback,1)];


            %Default reflection
            if isempty(obj.Reflect); obj.Reflect=0;end
        end





        function E = MoveField(obj,E,relative)
            E = MoveField@MBE.System.system(obj, E,relative);
            %Move backward fields in opposite direction
            E(end,obj.Nlines+1:end) = obj.Reflect(1) .* E(end,obj.ActiveBack);
            
        end
        
        function text = textIDsingle(obj,i)
            if i <= obj.Nlines
                text = textIDsingle@MBE.System.SM_Fwrd_Raman(obj,i);
            else
                order=obj.StokesList(i-obj.Nlines);
                if order == 0
                    text = 'Pump';
                elseif order>0
                    text = ['Back S ' num2str(order)];
                else
                    text = ['Back AS ' num2str(-order)];
                end
            end
        end

        function text = textQIDsingle(obj,i)
            BackQ=i-(obj.NQ-length(obj.StokesList));
            if BackQ <= 0
                text = textQIDsingle@MBE.System.SM_Fwrd_Raman(obj,i);
            else
                text = ['Q Back' num2str(obj.StokesList(BackQ))];
            end
        end

    end
end

