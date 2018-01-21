
function [Molar, Box, Solid] = constants(Calcite, Aragonite, Ca_fluid, Mg_fluid, C_fluid, Sr_fluid, Solid, alpha)

% Function that calculated constants used to initialize the carbonate 
% diagenesis model                         

global Fluid

% Molar Mass for different elements:
Molar.Mg = 24.305;                          % molar mass (g/mol)
Molar.Ca = 40.078;
Molar.C = 12.011;
Molar.O = 15.999;
Molar.Sr = 87.62;
 
% Sediment and box dimensions:
Box.vol = 100^3;                            % box volume (cm3/m3)
Box.Phi = 0.5;                              % Porosity (w. fraction)
Box.RhoS = 1.8;                             % Density of sediment (g/cm3)
Box.RhoF = 1.0125;                          % Density of fluid (g/cm3)
Box.sed = (1-Box.Phi)*Box.vol*Box.RhoS/1e3; % Sediment mass (kg)
Box.fluid = Box.Phi*Box.vol*Box.RhoF/1e3;   % Fluid mass (kg)

% Weight fraction of fluid in the system (cf. Banner and Hanson, 1990)
Box.F = (Box.Phi*Box.RhoF)/(Box.Phi*Box.RhoF + (1 - Box.Phi)*Box.RhoS);


%------------------- Initial composition of solid ------------------------%


% Elemental composition in calcite (CaCO3):

Solid.calc.init = Calcite*Box.sed;          % Initial mass of calcite (kg)

Solid.calc.Ca = 394500*Solid.calc.init...
    /(Molar.Ca*1e3);                        % Ca in calcite (mol)
Solid.calc.Mg = 5000*Solid.calc.init...
    /(Molar.Mg*1e3);                        % Mg in calcite (mol)
Solid.calc.C = 120000*Solid.calc.init...
    /(Molar.C*1e3);                         % C in calcite (mol)
Solid.calc.O = 480000*Solid.calc.init...
    /(Molar.O*1e3);                         % O in calcite (mol)
Solid.calc.Sr = 500*Solid.calc.init...
    /(Molar.Sr*1e3);

% Elemental composition in aragonite (CaCO3):

Solid.arag.init = Aragonite*Box.sed;        % Initial mass of aragonite (kg)
                                           
Solid.arag.Ca = 387500*Solid.arag.init...   % ppm converted to mol
    /(Molar.Ca*1e3); 
Solid.arag.Mg = 2500*Solid.arag.init...
    /(Molar.Mg*1e3);
Solid.arag.C = 120000*Solid.arag.init...
    /(Molar.C*1e3);
Solid.arag.O = 480000*Solid.arag.init...
    /(Molar.O*1e3);
Solid.arag.Sr = 10000*Solid.arag.init...
    /(Molar.Sr*1e3);

% For example: change Aragonite to High-Mg calcite: 

% High-Mg calcite:
% Solid.arag.Ca = 387000*Solid.arag.init/(Molar.Ca*1e3);  
% Solid.arag.Mg = 12000*Solid.arag.init/(Molar.Mg*1e3);
% Solid.arag.C = 120000*Solid.arag.init/(Molar.C*1e3);
% Solid.arag.O = 480000*Solid.arag.init/(Molar.O*1e3);
% Solid.arag.Sr = 1000*Solid.arag.init/(Molar.Sr*1e3);

% Elemental composition of dolomite in the solid (CaMg(CO3)2):

Dolomite = 0.000001;             
% Fraction of initial rock that is DIAGENETIC MINERAL (should be > 0)!

% Initial mass of diagenetic mineral (kg):
Solid.dolo.init = Dolomite*Box.sed;        

% Dolomite:
Solid.dolo.Ca = 217300*Solid.dolo.init/(Molar.Ca*1e3);
Solid.dolo.Mg = 131800*Solid.dolo.init/(Molar.Mg*1e3);
Solid.dolo.C = 130300*Solid.dolo.init/(Molar.C*1e3);
Solid.dolo.O = 520600*Solid.dolo.init/(Molar.O*1e3);
Solid.dolo.Sr = 50*Solid.dolo.init/(Molar.Sr*1e3);


%------------------- Isotopic composition of solid -----------------------%


% Isotopic composition of diagenetic mineral:
Solid.dolo.dCa = Solid.dolo.Ca*(Fluid.dCa + 1e3*log(alpha.Ca));    
Solid.dolo.dMg = Solid.dolo.Mg*(Fluid.dMg + 1e3*log(alpha.Mg));      
Solid.dolo.dC = Solid.dolo.C*(Fluid.dC + 1e3*log(alpha.C));    
Solid.dolo.dO = Solid.dolo.O*(Fluid.dO + 1e3*log(alpha.O)); 
Solid.dolo.dSr = Solid.dolo.Sr*(Fluid.dSr + 1e3*log(alpha.Sr));

% Solving d(M*delta)/dt requires initial conditions to be in units of
% mol*promil:

Solid.dCa = (Solid.calc.dCa*Solid.calc.Ca + Solid.arag.dCa*Solid.arag.Ca...
    + Solid.dolo.dCa);
Solid.dMg = (Solid.calc.dMg*Solid.calc.Mg + Solid.arag.dMg*Solid.arag.Mg...
    + Solid.dolo.dMg);
Solid.dC = (Solid.calc.dC*Solid.calc.C + Solid.arag.dC*Solid.arag.C...
    + Solid.dolo.dC);
Solid.dO = (Solid.calc.dO*Solid.calc.O + Solid.arag.dO*Solid.arag.O...
    + Solid.dolo.dO);
Solid.dSr = (Solid.calc.dSr*Solid.calc.Sr + Solid.arag.dSr*Solid.arag.Sr...
    + Solid.dolo.dSr);


%-------------------- Initial composition of fluid -----------------------%


% Fluid concentration:
Fluid.Ca = Ca_fluid*Box.fluid;              % Mass of Ca in the fluid (mol)
Fluid.Mg = Mg_fluid*Box.fluid;              % Mass of Mg in the fluid (mol) 
Fluid.C = C_fluid*Box.fluid;                % Mass of C (mol)
Fluid.O = 889e3*Box.fluid/(Molar.O*1e3);    % Mass of O (from ppm to mol)
Fluid.Sr = Sr_fluid*Box.fluid;              % Mass of Sr (mol) 

% Initial Mass*promil to solve ODE:
Fluid.dCa_init = Fluid.dCa*Fluid.Ca;       
Fluid.dMg_init = Fluid.dMg*Fluid.Mg;
Fluid.dC_init = Fluid.dC*Fluid.C;
Fluid.dO_init = Fluid.dO*Fluid.O;  
Fluid.dSr_init = Fluid.dSr*Fluid.Sr;

end