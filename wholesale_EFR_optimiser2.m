% load('ElecPrices.mat')

dt = 1; % hours
P_PV = 4; % MW
E_B = 4; % MWh
E_B_low = 0;
E_B_high = E_B;
P_B = 2; % MW
gridlim = 2; % MW
h1 = 7;
h2 = 19;
EFRrate = 250*P_B/8760; % MWh/h
footroom = EFRrate*(24-(h2-h1))+0.25*P_B;
R1=0.1;
eta = 0.95*((4-2*R1)/(4+2*R1)); % round-trip efficiency
etac = sqrt(eta);
etad = sqrt(eta);

load('n2exDATA.mat')
load('SgurrWindSolar.mat')

dtPV = 5/60; % hours
Pgen = zeros(size(priceblock15(:)));
p_elec = zeros(size(priceblock15(:)));
p_elec0 = priceblockTall(:);
offset = 365*24*2+1;

for nh = 1:24
    for nd = 1:365
        Pgen(nh+(nd-1)*24) = mean(Solar5Min( (dt/dtPV)*(24*(nd-1)+nh-1)+1 : (dt/dtPV)*(24*(nd-1)+nh) ));
        p_elec(nh+(nd-1)*24) = p_elec0(offset+ (24*(nd-1)+nh-1)+1 : offset+ (24*(nd-1)+nh) );
    end
end
% Pgen = P_PV*[Pgen(:,end) Pgen]/100; % MW
Pgen = P_PV*Pgen/100; % MW

% time series
PtoGrid = zeros(size(priceblock15(:)));
Pcharge = zeros(size(priceblock15(:)));
Estored = footroom*ones(size(priceblock15(:)));

nT = h2-h1;
pe = optimvar('pe',nT,'Type','integer','LowerBound',-gridlim,...
    'UpperBound',gridlim);
pc = optimvar('pc',nT,'LowerBound',0,'UpperBound',P_B);
pd = optimvar('pd',nT,'LowerBound',0,'UpperBound',P_B);
    
for nd = 1:365
    % 1 -  model EFR consuming battery SoC
    Estored_s_day = footroom*ones(1,24);
    Pcharge_s_day = zeros(1,24);
    PtoGrid_s_day = zeros(1,24);
    
    if h1>0
      if nd == 1
        Estored_s_day(1) = Estored_s_day(end)-EFRrate;
      else
        Estored_s_day(1) = Estored((nd-1)*24)-EFRrate;
      end
    end
    if h1>1
      for nh = 2:h1
        Estored_s_day(nh) = Estored_s_day(nh-1)-EFRrate;
      end
    elseif h1<1 % if day-mode begins at midnight
      E0 = Estored_s_day(1)+EFRrate;
    end
    
    % 2 - initiate charging/export
    if h1>0
      for nh = 1:h1
        % charge if sun rises before day-mode begins
        Pcharge_s_day(nh) = min([P_B,Pgen(nh+(nd-1)*24)]);
        Estored_s_day(nh) = Estored_s_day(nh) + Pcharge_s_day(nh)*dt*etac;
        if Estored_s_day(nh) > E_B % case: spillover charge
            Pcharge_s_day(nh) = Pcharge_s_day(nh) -(Estored_s_day(nh)-E_B)/(dt*etac);
            Estored_s_day(nh) = E_B;
        end
      end
      E0 = Estored_s_day(h1); % day-mode starting energy from previous hour
    end
    
    PGEN = Pgen(h1+(nd-1)*24+1:h2+(nd-1)*24);
    P_ELEC = p_elec(h1+(nd-1)*24+1:h2+(nd-1)*24);
    
    % 3 - optimise pe, pc, pd
    PowerBal = pc+pe-pd <= PGEN;
    
    Euplim = optimconstr(nT);
    Elowlim = optimconstr(nT);
    for t = 1:nT
        % must remain below upper charge limit
       Euplim(t) = E0 +sum(pc(1:t))*etac -sum(pd(1:t))/etad <= E_B_high;
        % must remain above lower charge limit
       Elowlim(t) = -(E0 +sum(pc(1:t))*etac -sum(pd(1:t))/etad) <= -E_B_low;
    end
    % in final hour, charge must be sufficient for night-mode operation
    Elowlim(nT) = -(E0 +sum(pc(1:nT))*etac -sum(pd(1:nT))/etad) <= -footroom;
    
    MaxProfit = sum(P_ELEC.*pe); % revenue to maximise
    
    dispatch = optimproblem('ObjectiveSense','maximize');
    dispatch.Objective = MaxProfit;
    dispatch.Constraints.PowerBal = PowerBal;
    dispatch.Constraints.Euplim = Euplim;
    dispatch.Constraints.Elowlim = Elowlim;

    options = optimoptions('intlinprog','Display','none','MaxTime',120);

    [dispatchsol,fval,exitflag,output] = solve(dispatch,options);
    
    % 4 - put optimised solutions into time series
    PtoGrid_s_day(h1+1:h2) = dispatchsol.pe;
    Pcharge_s_day(h1+1:h2) = dispatchsol.pc-dispatchsol.pd;
    for nh = h1+1:h2
        Estored_s_day(nh) = E0+etac*sum(dispatchsol.pc(1:nh-h1))...
                         -sum(dispatchsol.pd(1:nh-h1))/etad;
    end
        
    % 5 - work out charge/discharge for remainder of day
    if h2<24
      for nh = h2+1:24
        Estored_s_day(nh) = Estored_s_day(nh-1)-EFRrate;
        % charge if sun still hasn't set yet
        Pcharge_s_day(nh) = min([P_B,Pgen(nh+(nd-1)*24)]);
        Estored_s_day(nh) = Estored_s_day(nh) + Pcharge_s_day(nh)*dt*etac;
        if Estored_s_day(nh) > E_B % case: spillover charge
            Pcharge_s_day(nh) = Pcharge_s_day(nh) -(Estored_s_day(nh)-E_B)/(dt*etac);
            Estored_s_day(nh) = E_B;
        end
      end
    end
    
    Estored(1+(nd-1)*24:nd*24) = Estored_s_day;
    Pcharge(1+(nd-1)*24:nd*24) = Pcharge_s_day;
    PtoGrid(1+(nd-1)*24:nd*24) = PtoGrid_s_day;
    
end

t1 = 1/24:1/24:365;
Pgen1 = Pgen(:);
p_elec1 = p_elec(:);

figure
plot(t1,Pgen1,'r','linewidth',2)
hold on
grid on
plot(t1,p_elec1/10,'b')
plot(t1,PtoGrid(:),'.-','color',[0 0.6 0],'linewidth',1)
plot(t1,Pcharge(:),':','color',[0.7 0 0.6],'linewidth',1)
plot(t1,Estored(:),'color',[0.7 0.5 0])
legend('PV (MW)','price (10s £/MWh)','export (MW)','charge (MW)','stored (MWh)')
% xlim([159 165])

disp(['£',num2str((sum(sum(p_elec1.*PtoGrid-0.0235*abs(PtoGrid)))+7*(24-(h2-h1))*365*P_B)/40000),' mill.'])

% 
% disp(['average £',num2str(sum(sum(p_elec.*PtoGrid))/sum(sum(Pgen))),'/MWh'])
% disp(['£',num2str((sum(sum(p_elec.*PtoGrid))+7*(23-(h2-h1))*365*P_B)/40000),' mill.'])
% 
% disp('no battery')
% Pgen0 = Pgen;
% Pgen0(Pgen>gridlim) = gridlim;
% disp(['£',num2str(sum(sum(p_elec.*Pgen0))/40000),' mill.'])


