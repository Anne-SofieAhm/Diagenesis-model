
function [in, out] = fluxes6(calc, arag, dolo, fluid, flux, n)

% Function that solve all input and output fluxes to each box:

global alpha R M

% Define reaction rates for each mineral:
R_calc = R;                       
R_arag = R_calc;
R_dolo = R_calc*0.0;                    % re-dissolving diagenetic mineral

% Stoichiometry of the diagenetic mineral:
MgC = M;                                % for dolomite M = 0.5
CaC = (1 - MgC);                        % Ca flux is inverse of Mg
Kd = 0.025;                             % Distribution coeff. for Sr-Ca

%%
% Fluxes of precipitation of dolomite/diagenetic calcite:
% Stoichiometric scaled to dolomite CaMg(CO3)2

% Diagenetic Calcite dissolution/precipitation:
dia.calc.C = R_calc*calc.C;             % (mol/yr)
dia.calc.Mg = MgC*dia.calc.C;
dia.calc.Ca = CaC*dia.calc.C;
dia.calc.O = 3*dia.calc.C; 
dia.calc.Sr = R_calc*calc.Ca*Kd*fluid(n+1).Sr/fluid(n+1).Ca;

% Diagenetic Aragonite dissolution/precipitation:
dia.arag.C = R_arag*arag.C;             % (mol/yr)
dia.arag.Mg = MgC*dia.arag.C;
dia.arag.Ca = CaC*dia.arag.C;
dia.arag.O = 3*dia.arag.C; 
dia.arag.Sr = R_arag*Kd*arag.Ca*fluid(n+1).Sr/fluid(n+1).Ca;

% *Diagenetic dolomite dissolution/precipitation*
dia.dolo.C = R_dolo*dolo.C;             % (mol/yr)
dia.dolo.Mg = MgC*dia.dolo.C;
dia.dolo.Ca = CaC*dia.dolo.C;
dia.dolo.O = 3*dia.dolo.C;     
dia.dolo.Sr = R_dolo*Kd*dolo.Ca*fluid(n+1).Sr/fluid(n+1).Ca;


%% Sediment input and output fluxes (F_in - F_out)

% Fluxes of dissolution/precipitation is set by the rate constant (R) 

% Change in mass of primary calcite: 
in.calc.Ca = 0;                         % Precipitation of primary calcite
out.calc.Ca = R_calc*calc.Ca;           % Dissolution of primary calcite
in.calc.Mg = 0;                         % in mol/yr
out.calc.Mg = R_calc*calc.Mg;
in.calc.C = 0;
out.calc.C = R_calc*calc.C;
in.calc.O = 0;
out.calc.O = R_calc*calc.O;
in.calc.Sr = 0;
out.calc.Sr = R_calc*calc.Sr;

% Change in mass of primary aragonite:
in.arag.Ca = 0;                         % Precipitation of primary aragonite
out.arag.Ca = R_arag*arag.Ca;           % Dissolution of primary aragonite
in.arag.Mg = 0;                         % in mol/yr
out.arag.Mg = R_arag*arag.Mg;
in.arag.C = 0;
out.arag.C = R_arag*arag.C;
in.arag.O = 0;
out.arag.O = R_arag*arag.O;
in.arag.Sr = 0;
out.arag.Sr = R_arag*arag.Sr;

% Change in mass of diagenetic dolomite:
% Precipitation of the diagenetic mineral has to equal the dissolution of
% the solid in order to approximately conserve mass of the rock. 
% Hence, the dissolution flux of aragonite and calcite is equal to the 
% precipitation flux of dolomite, but scaled stoichimetric (by M).

in.dolo.Ca = dia.calc.Ca + dia.arag.Ca + dia.dolo.Ca; % [mol/yr]
out.dolo.Ca = R_dolo*dolo.Ca;
in.dolo.Mg = dia.calc.Mg + dia.arag.Mg + dia.dolo.Mg;
out.dolo.Mg = R_dolo*dolo.Mg;
in.dolo.C = dia.calc.C + dia.arag.C + dia.dolo.C;
out.dolo.C = R_dolo*dolo.C;
in.dolo.O = dia.calc.O + dia.arag.O + dia.dolo.O;
out.dolo.O = R_dolo*dolo.O;
in.dolo.Sr = dia.calc.Sr + dia.arag.Sr + dia.dolo.Sr;
out.dolo.Sr = R_dolo*dolo.Sr;


%% Fluid input and output fluxes:

% The components of the fluid is controled by initial seawater flushing 
% through the sediment, the continous dissolution of primary minerals, 
% and re- precipitation of diagenetic minerals. 

% fluid.Ca(n) refers to the composition of the previous box while
% fluid.Ca(n+1) referes to the updated fluid composition.

% Calcium in the fluid:
in.fluid.Ca = flux*fluid(n).Ca + R_calc*calc.Ca + R_arag*arag.Ca + R_dolo*dolo.Ca;
out.fluid.Ca = flux*fluid(n+1).Ca + dia.calc.Ca + dia.arag.Ca + dia.dolo.Ca;

% Magnesium in the fluid:
in.fluid.Mg = flux*fluid(n).Mg + R_calc*calc.Mg + R_arag*arag.Mg + R_dolo*dolo.Mg;
out.fluid.Mg = flux*fluid(n+1).Mg + dia.calc.Mg + dia.arag.Mg + dia.dolo.Mg;

% Carbon in the fluid:
in.fluid.C = flux*fluid(n).C + R_calc*calc.C + R_arag*arag.C + R_dolo*dolo.C;
out.fluid.C = flux*fluid(n+1).C + dia.calc.C + dia.arag.C + dia.dolo.C;

% Oxygen in the fluid:
in.fluid.O = flux*fluid(n).O + R_calc*calc.O + R_arag*arag.O + R_dolo*dolo.O;
out.fluid.O = flux*fluid(n+1).O + dia.calc.O + dia.arag.O + dia.dolo.O;

% Strontium in the fluid:
in.fluid.Sr = flux*fluid(n).Sr + R_calc*calc.Sr + R_arag*arag.Sr + R_dolo*dolo.Sr;
out.fluid.Sr = flux*fluid(n+1).Sr + dia.calc.Sr + dia.arag.Sr + dia.dolo.Sr;

%% Isotopic composition

% Mass [mol] of the individual isotopes in the solid is not constant so
% dM*delta/dt has to be solved without isolating delta! This means that 
% all delta values that are solve for has units of mol*delta and has to be
% dived by mol to get delta. This is accounted for in all equations so that
% the units of flux always is mol/yr.

% Set fractionation factors:
dia.d.Ca = fluid(n+1).dCa + 1e3*log(alpha.Ca);    
dia.d.Mg = fluid(n+1).dMg + 1e3*log(alpha.Mg);
dia.d.C = fluid(n+1).dC + 1e3*log(alpha.C);
dia.d.O = fluid(n+1).dO + 1e3*log(alpha.O);
dia.d.Sr = fluid(n+1).dSr + 1e3*log(alpha.Sr);

%%
% Isotopes of the dolomite component
% Only changing in the case where the diagenetic mineral / dolomite is
% being re-dissolved by setting R_dolo (see top)

in.dolo.dCa = (dia.calc.Ca + dia.arag.Ca + dia.dolo.Ca)*dia.d.Ca;
out.dolo.dCa = R_dolo*dolo.Ca*dolo.dCa;

in.dolo.dMg = (dia.calc.Mg + dia.arag.Mg + dia.dolo.Mg)*dia.d.Mg;
out.dolo.dMg = R_dolo*dolo.Mg*dolo.dMg;

in.dolo.dC = (dia.calc.C + dia.arag.C + dia.dolo.C)*dia.d.C;
out.dolo.dC = R_dolo*dolo.C*dolo.dC;

in.dolo.dO = (dia.calc.O + dia.arag.O + dia.dolo.O)*dia.d.O;
out.dolo.dO = R_dolo*dolo.O*dolo.dO;

in.dolo.dSr = (dia.calc.Sr + dia.arag.Sr + dia.dolo.Sr)*dia.d.Sr;
out.dolo.dSr = R_dolo*dolo.Sr*dolo.dSr;

%% Ca isotopes:

% _Ca isotopes in the solid (mol*promil):_
in.solid.dCa = (dia.calc.Ca + dia.arag.Ca + dia.dolo.Ca)*dia.d.Ca;
out.solid.dCa = R_calc*calc.Ca*calc.dCa + R_arag*arag.Ca*arag.dCa + R_dolo*dolo.Ca*dolo.dCa;

% _Ca isotopes in the fluid (mol*promil):_
in.fluid.dCa = flux*fluid(n).Ca*fluid(n).dCa + R_calc*calc.Ca*calc.dCa +...
    R_arag*arag.Ca*arag.dCa + R_dolo*dolo.Ca*dolo.dCa;
out.fluid.dCa = flux*fluid(n+1).Ca*fluid(n+1).dCa +...
    (dia.calc.Ca + dia.arag.Ca + dia.dolo.Ca)*dia.d.Ca;

%% Mg isotopes:

% _Mg isotopes in the solid:_
in.solid.dMg = (dia.calc.Mg + dia.arag.Mg + dia.dolo.Mg)*dia.d.Mg;
out.solid.dMg = R_calc*calc.Mg*calc.dMg + R_arag*arag.Mg*arag.dMg + R_dolo*dolo.Mg*dolo.dMg;
 
% _Mg isotopes in the fluid:_
in.fluid.dMg = flux*fluid(n).Mg*fluid(n).dMg + R_calc*calc.Mg*calc.dMg +...
    R_arag*arag.Mg*arag.dMg + R_dolo*dolo.Mg*dolo.dMg;
out.fluid.dMg = flux*fluid(n+1).Mg*fluid(n+1).dMg +...
    (dia.calc.Mg + dia.arag.Mg + dia.dolo.Mg)*dia.d.Mg;

%% C isotopes:

% _C isotopes in the solid:_
in.solid.dC = (dia.calc.C + dia.arag.C + dia.dolo.C)*dia.d.C;
out.solid.dC = R_calc*calc.C*calc.dC + R_arag*arag.C*arag.dC + R_dolo*dolo.C*dolo.dC;

% _C isotopes in the fluid:_
in.fluid.dC = flux*fluid(n).C*fluid(n).dC + R_calc*calc.C*calc.dC +...
    R_arag*arag.C*arag.dC + R_dolo*dolo.C*dolo.dC;
out.fluid.dC = flux*fluid(n+1).C*fluid(n+1).dC +...
    (dia.calc.C + dia.arag.C + dia.dolo.C)*dia.d.C;

%% O isotopes:

% _O isotopes in the solid:_
in.solid.dO = (dia.calc.O + dia.arag.O + dia.dolo.O)*dia.d.O;
out.solid.dO = R_calc*calc.O*calc.dO + R_arag*arag.O*arag.dO + R_dolo*dolo.O*dolo.dO;

% _O isotopes in the fluid:_
in.fluid.dO = flux*fluid(n).O*fluid(n).dO + R_calc*calc.O*calc.dO +...
    R_arag*arag.O*arag.dO + R_dolo*dolo.O*dolo.dO;
out.fluid.dO = flux*fluid(n+1).O*fluid(n+1).dO +...
    (dia.calc.O + dia.arag.O + dia.dolo.O)*dia.d.O;

%% Sr isotopes:

% _Sr isotopes in the solid:_
in.solid.dSr = (dia.calc.Sr + dia.arag.Sr + dia.dolo.Sr)*dia.d.Sr;
out.solid.dSr = R_calc*calc.Sr*calc.dSr + R_arag*arag.Sr*arag.dSr + R_dolo*dolo.Sr*dolo.dSr;

% _Sr isotopes in the fluid(n+1):_
in.fluid.dSr = flux*fluid(n).Sr*fluid(n).dSr + R_calc*calc.Sr*calc.dSr +...
    R_arag*arag.Sr*arag.dSr + R_dolo*dolo.Sr*dolo.dSr;
out.fluid.dSr = flux*fluid(n+1).Sr*fluid(n+1).dSr +...
    (dia.calc.Sr + dia.arag.Sr + dia.dolo.Sr)*dia.d.Sr;


end