function dy = dJ6dt(t,J,m,v)
global Fluid

box = length(m);
dy = zeros(box*v,1);

%% box 1

% Initial conditions:
fluid(1).Ca = Fluid.Ca;     % Initial input in box 1 is always sea-water (mol)
fluid(1).Mg = Fluid.Mg;
fluid(1).C = Fluid.C;
fluid(1).O = Fluid.O;
fluid(1).Sr = Fluid.Sr;
fluid(1).dCa = Fluid.dCa;
fluid(1).dMg = Fluid.dMg;
fluid(1).dC = Fluid.dC;
fluid(1).dO = Fluid.dO;
fluid(1).dSr = Fluid.dSr;

%% rest of boxes
for n = 1:box
    
    i = v*n - v;                            % vector with ODE results
                                            % MatLab refuses to work with
                                            % matrixes in an ODE ?!?
    
% Mineral mixing / Mass of solid:
calc.Ca = J(i+1);                           % Ca in calcite (mol)
calc.Mg = J(i+2);                           % Mg in calcite (mol)
calc.C = J(i+3);                            % C in calcite (mol)
calc.O = J(i+4);                            % O in calcite (mol)
calc.Sr = J(i+5);                           % Sr in calcite (mol)

arag.Ca = J(i+6);                           % Ca in aragonite (mol)                  
arag.Mg = J(i+7);                           % Mg in aragonite (mol)
arag.C = J(i+8);                            % C in aragonite (mol)
arag.O = J(i+9);                            % O in aragonite (mol)
arag.Sr = J(i+10);                          % Sr in aragonite (mol)

dolo.Ca = J(i+11);                          % Ca in dolomite (mol)
dolo.Mg = J(i+12);                          % Mg in dolomite (mol)
dolo.C = J(i+13);                           % C in dolomite (mol)
dolo.O = J(i+14);                           % O in dolomite (mol)
dolo.Sr = J(i+15);                          % Sr in dolomite (mol)

% Isotope composition of solid:
calc.dCa = J(i+16);                         % d44Ca in calcite (per mil)
calc.dMg = J(i+17);                         % d26Mg in calcite (per mil)
calc.dC = J(i+18);                          % d13C in calcite (per mil)
calc.dO = J(i+19);                          % d18O in calcite (per mil)
calc.dSr = J(i+20);                         % 86Sr/87Sr in calcite

arag.dCa = J(i+21);                         % d44Ca in aragonite (per mil)
arag.dMg = J(i+22);                         % d26Mg in aragonite (per mil)
arag.dC = J(i+23);                          % d13C in aragonite (per mil)
arag.dO = J(i+24);                          % d18O in aragonite (per mil)
arag.dSr = J(i+25);                         % 86Sr/87Sr in aragonite

dolo.dCa = J(i+26)/dolo.Ca;                 % d44Ca in dolomite (per mil)
dolo.dMg = J(i+27)/dolo.Mg;                 % d26Mg in dolomite (per mil)
dolo.dC = J(i+28)/dolo.C;                   % d13C in dolomite (per mil)
dolo.dO = J(i+29)/dolo.O;                   % d18O in dolomite (per mil)
dolo.dSr = J(i+30)/dolo.Sr;                 % 86Sr/87Sr in dolomite

% Mass of the fluid:
fluid(n+1).Ca = J(i+36);                    % Ca in fluid (mol)
fluid(n+1).Mg = J(i+37);                    % Mg in fluid (mol)
fluid(n+1).C = J(i+38);                     % C in fluid (mol)
fluid(n+1).O = J(i+39);                     % O in fluid (mol)
fluid(n+1).Sr = J(i+40);                    % Sr in fluid (mol)

% Isotope composition of fluid:
fluid(n+1).dCa = J(i+41)/fluid(n+1).Ca;     % d44Ca in fluid (mol*promil/mol)
fluid(n+1).dMg = J(i+42)/fluid(n+1).Mg;     % d26Mg in fluid (mol*promil/mol)
fluid(n+1).dC = J(i+43)/fluid(n+1).C;       % d13C in fluid (mol*promil/mol)
fluid(n+1).dO = J(i+44)/fluid(n+1).O;       % d18O in fluid (mol*promil/mol)
fluid(n+1).dSr = J(i+45)/fluid(n+1).Sr;     % 86Sr/87Sr in fluid (mol*promil/mol)

% Flux
flux = J(i+46);                             % Fluid flow rate

%%
% Call function that solve all input and output fluxes to each box:

 [in, out] = fluxes6(calc, arag, dolo, fluid, flux, n);

%%
% Solve ODE:
% (dM/dt = input fluxes - output fluxes)

dy(i+1) = in.calc.Ca - out.calc.Ca;
dy(i+2) = in.calc.Mg - out.calc.Mg;
dy(i+3) = in.calc.C - out.calc.C;
dy(i+4) = in.calc.O - out.calc.O;
dy(i+5) = in.calc.Sr - out.calc.Sr;

dy(i+6) = in.arag.Ca - out.arag.Ca;
dy(i+7) = in.arag.Mg - out.arag.Mg;
dy(i+8) = in.arag.C - out.arag.C;
dy(i+9) = in.arag.O - out.arag.O;
dy(i+10) = in.arag.Sr - out.arag.Sr;

dy(i+11) = in.dolo.Ca - out.dolo.Ca;
dy(i+12) = in.dolo.Mg - out.dolo.Mg;
dy(i+13) = in.dolo.C - out.dolo.C;
dy(i+14) = in.dolo.O - out.dolo.O;
dy(i+15) = in.dolo.Sr - out.dolo.Sr;

dy(i+16) = 0;     % Constant isotopic end-members of the initial composition
dy(i+17) = 0;
dy(i+18) = 0;
dy(i+19) = 0;
dy(i+20) = 0;
dy(i+21) = 0;
dy(i+22) = 0;
dy(i+23) = 0;
dy(i+24) = 0;
dy(i+25) = 0;

dy(i+26) = in.dolo.dCa - out.dolo.dCa;
dy(i+27) = in.dolo.dMg - out.dolo.dMg;
dy(i+28) = in.dolo.dC - out.dolo.dC;
dy(i+29) = in.dolo.dO - out.dolo.dO;
dy(i+30) = in.dolo.dSr - out.dolo.dSr;

dy(i+31) = (in.solid.dCa - out.solid.dCa); % NB! notice not divided by mass
dy(i+32) = in.solid.dMg - out.solid.dMg;
dy(i+33) = in.solid.dC - out.solid.dC;
dy(i+34) = in.solid.dO - out.solid.dO;
dy(i+35) = in.solid.dSr - out.solid.dSr;

dy(i+36) = in.fluid.Ca - out.fluid.Ca;
dy(i+37) = in.fluid.Mg - out.fluid.Mg;
dy(i+38) = in.fluid.C - out.fluid.C;
dy(i+39) = in.fluid.O - out.fluid.O;
dy(i+40) = in.fluid.Sr - out.fluid.Sr;

dy(i+41) = (in.fluid.dCa - out.fluid.dCa);
dy(i+42) = in.fluid.dMg - out.fluid.dMg;
dy(i+43) = in.fluid.dC - out.fluid.dC;
dy(i+44) = in.fluid.dO - out.fluid.dO;
dy(i+45) = in.fluid.dSr - out.fluid.dSr;

dy(i+46) = 0;     % Fluid flow rate kept constant

end

%t

end
