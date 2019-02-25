%%                  Carbonate dissolution-precipitation model
%%
% V6.0: Dissolution of carbonate and precipitation of diagenetic mineral
% 08-04-2015 by Anne-Sofie C. Ahm 
% Revised 11-30-2015 to include Sr
% Published in Ahm et al. (2018) GCA Vol. XX No. XX pp. 

clear all; clc;
global Fluid alpha R M

%%                          Model Variables:
%%

%------------------- Composition of primary rock:-------------------------%


% Initial mineralogy:
Calcite = 0.01;              % Fraction of calcite (has to be > 0)
Aragonite = 0.99;            % Fraction of aragonite (has to be > 0)

% Calcite:                   Aragonite:              
Solid.calc.dCa = -1.1;       Solid.arag.dCa = -1.6;     % d44Ca value      
Solid.calc.dMg = -3.8;       Solid.arag.dMg = -3.8;     % d26Mg value 
Solid.calc.dC = 3.2;         Solid.arag.dC = 5;         % d13C value
Solid.calc.dO = -2;          Solid.arag.dO = 0;         % d18O value
Solid.calc.dSr = 0.70900;    Solid.arag.dSr = 0.70900;  % 87Sr/86Sr ratio


%----------------- Composition of diagenetic fluid:-----------------------%


% Concentration:
Ca_fluid = 0.01028;          % Ca in fluid (mol kg-1)
Mg_fluid = 0.05282;          % Mg in fluid (mol kg-1)
C_fluid = 0.002;             % C in fluid (mol kg-1)
Sr_fluid = 0.00009;          % Sr in fluid (mol kg-1)

% Isotopic composition:
Fluid.dCa = 0;               % d44Ca_fluid value
Fluid.dMg = -0.8;            % d26Mg_fluid value  
Fluid.dC = -2;               % d13C_fluid value   
Fluid.dO = -29.5;            % d18O_fluid value (normalized to V-PDB)
Fluid.dSr = 0.709250;        % 86Sr/87Sr_fluid value

% Fractionation factors of diagenetic mineral: 
alpha.Ca = 1.0;              % d44Ca_alpha 
alpha.O = 1.0289;            % d18O_alpha 
alpha.C = 1.001;             % d13C_alpha
alpha.Mg = 0.998;            % d26Mg_alpha 
alpha.Sr = 1;                % 87Sr/86Sr_alpha (should be 1)   


%--------------------- Diagenetic Lenght Scale ---------------------------%


u = 0.10;                    % Advection Rate (m yr-1)                
R = 10e-6;                   % Reaction Rate (yr-1)
time = 1e7;                  % Timespan modeled (yrs)
box = 20;                    % Number of boxes (i.e sediment depth) (m)

                             % NB! Dependinding on the reaction rate the 
                             % problem becomes too stiff for the ODE solver 
                             % if there are too many boxes relative to the 
                             % advection/reaction rate

%------------------ Stoichiometry of diagenetic mineral ------------------%
M =  0.5;                    
                             % Mg/Ca ratio of diagenetic mineral 
                             % (0.5 = dolomite, 0.035 = low-Mg calcite)                   

                             
%%                           Constants:
%%

% Function that calls constants:
        
[Molar, Box, Solid] = constants(Calcite, Aragonite, Ca_fluid, Mg_fluid,...
    C_fluid, Sr_fluid, Solid, alpha);

                            % Molar = elemental molar mass
                            % Box = box volume and porosity
                            % Solid = Initial composition of bulk rock
                            % (e.g. concentration of elements in primary
                            % minerals)
                             
                            
%%                  Model Input Vector
%%   

% Initial value vector (to be passed to ODE):

J = [Solid.calc.Ca; Solid.calc.Mg; Solid.calc.C; Solid.calc.O;...
    Solid.calc.Sr; Solid.arag.Ca; Solid.arag.Mg; Solid.arag.C;...
    Solid.arag.O; Solid.arag.Sr; Solid.dolo.Ca; Solid.dolo.Mg;...
    Solid.dolo.C; Solid.dolo.O; Solid.dolo.Sr; Solid.calc.dCa;...
    Solid.calc.dMg; Solid.calc.dC; Solid.calc.dO; Solid.calc.dSr;...
    Solid.arag.dCa; Solid.arag.dMg; Solid.arag.dC; Solid.arag.dO;...
    Solid.arag.dSr; Solid.dolo.dCa; Solid.dolo.dMg; Solid.dolo.dC;...
    Solid.dolo.dO; Solid.dolo.dSr; Solid.dCa; Solid.dMg; Solid.dC;...
    Solid.dO; Solid.dSr; Fluid.Ca; Fluid.Mg; Fluid.C; Fluid.O;...
    Fluid.Sr; Fluid.dCa_init; Fluid.dMg_init; Fluid.dC_init;...
    Fluid.dO_init; Fluid.dSr_init; u];

 
%%                  Solve finite difference scheme
%%

% Set initial boundary conditions:
tspan = [0 time];            % time span (yr)
m = 1:box;                   % vector with number of boxes

for i = 1:length(m)          % Loop passing initial conditions to each box
        
    Jt0(:,i) = J(:,1); 

end

v = length(Jt0(:,1));        % Number of variables in initial vector

% Solve for change over time (dM/dt) by calling ODE solver (dJ6dt):
tic, [t,J] = ode15s(@(t,J)dJ6dt(t,J,m,v),tspan,Jt0); toc           
                                                      
% Loop to compress all output data into a single 3d matrix (JPrime) that
% can be saved ('save JPrime'):
for i = 1:length(m) 

        j = v*i - v;
        a = j + 1; c = j + v;                                
        JPrime(:,:,i) = J(1:length(t),a:c)';  % Model output 3D-matrix        
        
end


%%                           Model output
%%

% Extract variables from model output 3d matrix (JPrime)

for i = 1:length(t)
    for j = 1:length(m)

% Elemental mass in calcite (mol):
Solid.calc.Ca(i,j) = JPrime(1,i,j);          
Solid.calc.Mg(i,j) = JPrime(2,i,j);          
Solid.calc.C(i,j) = JPrime(3,i,j);            
Solid.calc.O(i,j) = JPrime(4,i,j);           
Solid.calc.Sr(i,j) = JPrime(5,i,j);

% Elemental mass in aragonite (mol):
Solid.arag.Ca(i,j) = JPrime(6,i,j);
Solid.arag.Mg(i,j) = JPrime(7,i,j);
Solid.arag.C(i,j) = JPrime(8,i,j);
Solid.arag.O(i,j) = JPrime(9,i,j);
Solid.arag.Sr(i,j) = JPrime(10,i,j);

% Elemental mass in the diagenetic mineral (mol):
Solid.dolo.Ca(i,j) = JPrime(11,i,j);
Solid.dolo.Mg(i,j) = JPrime(12,i,j);
Solid.dolo.C(i,j) = JPrime(13,i,j);
Solid.dolo.O(i,j) = JPrime(14,i,j);
Solid.dolo.Sr(i,j) = JPrime(15,i,j);

% Total elemental mass (mol) in the bulk sediment: 
Total_Ca(i,j) =...
    (Solid.calc.Ca(i,j) + Solid.arag.Ca(i,j) + Solid.dolo.Ca(i,j));
Total_Mg(i,j) =...
    (Solid.calc.Mg(i,j) + Solid.arag.Mg(i,j) + Solid.dolo.Mg(i,j));
Total_C(i,j) =...
    (Solid.calc.C(i,j) + Solid.arag.C(i,j) + Solid.dolo.C(i,j));
Total_O(i,j) =...
    (Solid.calc.O(i,j) + Solid.arag.O(i,j) + Solid.dolo.O(i,j));
Total_Sr(i,j) =...
    (Solid.calc.Sr(i,j) + Solid.arag.Sr(i,j) + Solid.dolo.Sr(i,j));

% Calculate elemetal ratios in the bulk sediment (mmol/mol):
Solid.Sr_Ca = Total_Sr./(Total_Ca)*1e3;
Solid.Mg_Ca = Total_Mg./Total_Ca*1e3;

% Total mass (mol) of the bulk sediment:
FWs(i,j) = Total_Ca(i,j) + Total_Mg(i,j) + Total_C(i,j)...
    + Total_O(i,j) + Total_Sr(i,j);

% Isotopic composition of the diagenetic mineral:
Solid.dolo.dCa(i,j) = JPrime(26,i,j)/Solid.dolo.Ca(i,j);
Solid.dolo.dMg(i,j) = JPrime(27,i,j)/Solid.dolo.Mg(i,j);
Solid.dolo.dC(i,j) = JPrime(28,i,j)/Solid.dolo.C(i,j);
Solid.dolo.dO(i,j) = JPrime(29,i,j)/Solid.dolo.O(i,j);
Solid.dolo.dSr(i,j) = JPrime(30,i,j)/Solid.dolo.Sr(i,j);

% Isotopic composition of the bulk sediment:
Solid.dCa(i,j) = JPrime(31,i,j)/Total_Ca(i,j);               
Solid.dMg(i,j) = JPrime(32,i,j)/Total_Mg(i,j);
Solid.dC(i,j) = JPrime(33,i,j)/Total_C(i,j);
Solid.dO(i,j) = JPrime(34,i,j)/Total_O(i,j);
Solid.dSr(i,j) = JPrime(35,i,j)/Total_Sr(i,j);

% Elemental mass (mol) of the diagenetic fluid:
Fluid.Ca(i,j) = JPrime(36,i,j);
Fluid.Mg(i,j) = JPrime(37,i,j);
Fluid.C(i,j) = JPrime(38,i,j);
Fluid.O(i,j) = JPrime(39,i,j);
Fluid.Sr(i,j) = JPrime(40,i,j);

% Total mass (mol) of the diagenetic fluid:
FWf(i,j) = Fluid.Ca(i,j) + Fluid.Mg(i,j) + Fluid.C(i,j)...
    + Fluid.O(i,j) + Fluid.Sr(i,j);

% Isotopic composition of the diagenetic fluid:
Fluid.dCa(i,j) = JPrime(41,i,j)/Fluid.Ca(i,j);
Fluid.dMg(i,j) = JPrime(42,i,j)/Fluid.Mg(i,j);
Fluid.dC(i,j) = JPrime(43,i,j)/Fluid.C(i,j);
Fluid.dO(i,j) = JPrime(44,i,j)/Fluid.O(i,j);
Fluid.dSr(i,j) = JPrime(45,i,j)/Fluid.Sr(i,j);

% Cummulative fluid-to-rock ratio:
WR(i) = (t(i)*u*Box.fluid + Box.fluid)/Box.sed; 

% Procent alteration (exponential decay):
Solid.procent(i,j) = exp(-R*t(i));

    end 
end
    

%%                           Figures
%%


% -- Change in bulk rock chemsitry with increasing fluid-to-rock ratios --%


figure

suptitle('Diagenesis over time')

subplot(7,1,1)

semilogx(WR,Solid.procent(:,end)); hold on
semilogx(WR,1-Solid.procent(:,end));
    xlim([1e1 WR(end)]); ylabel('%')
    
subplot(7,1,2)

semilogx(WR,Solid.dC(:,[1,end])); hold on
    xlim([1e1 WR(end)]); ylabel('\delta^1^3C')

subplot(7,1,3)

semilogx(WR,Solid.dO(:,[1,end])); hold on
    xlim([1e1 WR(end)]); ylabel('\delta^1^8O')

subplot(7,1,4)

semilogx(WR,Solid.dCa(:,[1,end])); hold on
    xlim([1e1 WR(end)]); ylabel('\delta^4^4Ca')

subplot(7,1,5)

semilogx(WR,Solid.dMg(:,[1,end])); hold on
    xlim([1e1 WR(end)]); ylabel('\delta^2^6Mg')

subplot(7,1,6)

semilogx(WR,Solid.Sr_Ca(:,[1,end])); hold on
    xlim([1e1 WR(end)]); ylabel('Sr/Ca')

subplot(7,1,7)

semilogx(WR,Solid.dSr(:,[1,end])); hold on
    xlim([1e1 WR(end)]); ylabel('^{87}Sr/^{86}Sr')
    xlabel('Cumulative fluid-to-rock ratio')
     legend('box 1','box n')

   
%--- Model phase space defined by covariation between pairs of proxies ---%


figure

suptitle('Model phase space')

colr = [0.7 0.7 0.7]; 

subplot(2,2,1)

plot(Solid.dCa(:,[1,end]),Solid.Sr_Ca(:,[1,end]),'Color','k'); hold on
plot(Solid.dCa(end,:),Solid.Sr_Ca(end,:),'Color',[0.1 0.7 1],'LineWidth',1);
    contour(Solid.dCa,Solid.Sr_Ca,Solid.procent,'-','Color',[0.8 0.8 0.8]);
    xlabel('\delta^4^4Ca'); 
    ylabel('Sr/Ca (mmol/mol)'); 

subplot(2,2,2)

plot(Solid.dCa(:,[1,end]),Solid.dO(:,[1,end]),'Color','k'); hold on
plot(Solid.dCa(end,:),Solid.dO(end,:),'Color',[0.1 0.7 1],'LineWidth',1);
    contour(Solid.dCa,Solid.dO,Solid.procent,'-','Color',[0.8 0.8 0.8]);
    xlabel('\delta^4^4^/^4^0Ca'); 
    ylabel('\delta^1^8O'); 

subplot(2,2,3)

plot(Solid.dCa(:,[1,end]),Solid.dC(:,[1,end]),'Color','k'); hold on
plot(Solid.dCa(end,:),Solid.dC(end,:),'Color',[0.1 0.7 1],'LineWidth',1);
    contour(Solid.dCa,Solid.dC,Solid.procent,'-','Color',[0.8 0.8 0.8]);
    xlabel('\delta^4^4Ca');
    ylabel('\delta^1^3C');

subplot(2,2,4)

plot(Solid.dCa(:,[1,end]),Solid.dMg(:,[1,end]),'Color','k'); hold on
plot(Solid.dCa(end,:),Solid.dMg(end,:),'Color',[0.1 0.7 1],'LineWidth',1);
    contour(Solid.dCa,Solid.dMg,Solid.procent,'-','Color',[0.8 0.8 0.8]);
    xlabel('\delta^4^4Ca'); 
    ylabel('\delta^2^6Mg'); 



