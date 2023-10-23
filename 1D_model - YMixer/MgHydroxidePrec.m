%% Calculate population balance
clear;
close all;


%% Saving settings section
figureName    =  "figName";
isSavePng     = 0;                 % 1 saves, 0 doesn't save
%% - - - - - - - - - - - - - - - - CaseSetup - - - - - - - - - - - - - - - - %%
%% - - - - - - - - - - - - - -Structure definition- - - - - - - - - - - - - -%%
conc = 1;

s.kv=pi/6;
s.k_sp=power(10,-10.88);
s.rhoCry=2344.6;                               % kg/m^3
s.PMCry= 58.32;                                % g/mol

s.M_MgCl2=conc;                             % concentration @ inlet (mixtureFrac=0) M
s.M_NaOH=2*s.M_MgCl2;                          % concentration @ inlet (mixtureFrac=1) M
s.alpha=0.5;
s.nu=1e-6;                                     % m^2/s

s.nucleateSize=1e-9;                           % critical length, m
s.N=3;                                         % number of nodes
s.activity=true;

s.A_=1.49e26;                                  % (particles number)/(m^3*s)
s.B_=301;                                      % /

s.A_het=7.4e14;                                % (particles number)/(m^3*s)
s.B_het=30;                                    % /

s.kg_=10^(-9.6);                               % m/s
s.g_=1.001;                                    % /

s.diffContr=false;                             % enable diffusion controlled growth

s.kb=1.38e-23;                                 % J/K
s.T=(20+273.15);                               % K
s.mu=1e-3;                                     % Pa*s
s.c_hd=1/8*power((8*pi)/15,0.5);               % /
s.C_adjust=0.858;                              % C adjusted for hydro-dynamic
s.C_adjust2=s.C_adjust;                        % C adjusted for brownian

s.aggregation=true;                            % enable aggregation kernels
s.aggregEfficiency=true;
s.Ap=5.897;                                    % N/m^2 (aggregation efficient)

s.Cphi=2;                                      % microMixing parameter
s.variance0=0.25;                              % /
s.velocity=17.73;                              % m/s (fatto per 0.04 m)

% - - - - - - - - - - - - - -Structure definition- - - - - - - - - - - - - -%
% - - - - - - - - - - - - - - - - CaseSetup - - - - - - - - - - - - - - - - %

%% Call ODE for microMixing
options     =odeset('AbsTol',1e-10);
[t_]=microMixing(s);
[s.tVar,s.variance]=ode45(@varEvolution,[t_(1) t_(end)],s.variance0,options,s);
s.tend=s.tVar(end);                   % sec

%% Call for PBM solution
[tSol,ySol]=PBM(s);

%% Call for plotting solution
sizes = MgOH2Plot(tSol,ySol,s).*1e9 % diameter d10, 21, 32, 43 (sizes are in nanometers)