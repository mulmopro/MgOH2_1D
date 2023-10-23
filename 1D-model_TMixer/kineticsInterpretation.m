%% kinetics interpretation off-line
clear;
close all;
clc;
global nucleateSize velocity

% Primary nucleation, molecular growth, aggregation, expAggrEff
figureName    =  "primaryNucleation";

% Homo. nucl.
A=10^(26.1593);
B=301.44;

% Het. nucl.
Ahet=10^(14.8173);
Bhet=30.35;

S=linspace(10,1e6,1e7);

J=A.*exp(-B./power(log(S+1),2));
Jhet=Ahet.*exp(-Bhet./power(log(S+1),2));

% figure('DefaultAxesFontSize',16)
% semilogy(1./power(log(S+1),2),J,'LineWidth',2), hold on, semilogy(1./power(log(S+1),2),Jhet,'LineWidth',2)
% xlabel("$\frac{1}{ln^{2}(S+1)}$",'Interpreter','latex','FontSize',20)
% ylabel("ln(J)",'Interpreter','latex','FontSize',20)
% set(gcf,'color','w');

figure('DefaultAxesFontSize',10)
loglog(S,J,'k','LineWidth',2), hold on, loglog(S,Jhet,'--k','LineWidth',2)
xlabel("Supersaturation",'Interpreter','latex','FontSize',12)
ylabel("Primary nucleation rate",'Interpreter','latex','FontSize',12)
xticks([1,10,100,1000,10000,100000,1000000])
xlim([10,1000000])
set(gcf,'color','w');
set(gcf,'units','inches','position',[1,1,3.25,2.5])

%% Molecular growth
S=linspace(10,1e7,1e5);

% molecular growth
kg=10^(-9.6);                              % m/s
g=1.001;                                   % /

G=kg.*S.^g;

nucleateSize=1e-9; velocity=12; k_sp_=10^(-10.88); gammapm=0.6; PMCry_=58.34; rhoCry_=2344; T_=25+273; t=1e-3; nu_=1e-6;
L_=[1.00001e-9 ;1e-8; 1e-7];

G_DC=zeros(length(L_),length(S));
for j=1:length(L_)
    L=L_(j)
for i=1:length(S)
    if S(i)>0 && L>nucleateSize
        D=1e-9;                            % m^2/s 
        interTension=0.095;                % J/m^2
        R=8.314;                           % J/(mol*K)
        C_sInf=(k_sp_/(4*gammapm))^(1/3);
        v_mol=PMCry_/rhoCry_*0.001;        %m^3/mol
        expFactor=2*interTension*v_mol./(R*T_*L);
        C_sCurv=C_sInf*exp(expFactor);
        Re=epsilon(t).^(1/3)*L.^(4/3)./nu_;
        Sc=nu_./D;
        Sh=2+0.52*Re.^0.52*Sc.^(1/3);
        
        G_DC(j,i)=2*(Sh*D./L).*v_mol.*C_sCurv*S(i);
    else
        G_DC(j,i)=0;
    end
end
end

figure('DefaultAxesFontSize',16)
loglog(S,G,'LineWidth',2), hold on, loglog(S,G_DC,'LineWidth',2)
xlabel("S",'Interpreter','latex','FontSize',20)
ylabel("G, $\frac{m}{s}$",'Interpreter','latex','FontSize',20)
legend('Power Law','Diffusion controlled, L_{0}=1 nm','Diffusion controlled, L_{0}=10 nm','Diffusion controlled, L_{0}=100 nm')
set(gcf,'color','w');