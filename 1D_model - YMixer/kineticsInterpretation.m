%% kinetics interpretation off-line
clear;
close all;
clc;

% Primary nucleation, molecular growth, aggregation, expAggrEff

% Homo. nucl.
A=10^(26.1593);
B=301.44;

% Het. nucl.
Ahet=10^(14.8173);
Bhet=30.35;

S=linspace(1,1e7,1e7);

J=A.*exp(-B./power(log(S+1),2));
Jhet=Ahet.*exp(-Bhet./power(log(S+1),2));

% figure('DefaultAxesFontSize',16)
% semilogy(1./power(log(S+1),2),J,'LineWidth',2), hold on, semilogy(1./power(log(S+1),2),Jhet,'LineWidth',2)
% xlabel("$\frac{1}{ln^{2}(S+1)}$",'Interpreter','latex','FontSize',20)
% ylabel("ln(J)",'Interpreter','latex','FontSize',20)
% set(gcf,'color','w');

figure('DefaultAxesFontSize',10)
loglog(S,J,'k-','LineWidth',2), hold on, loglog(S,Jhet,'k--','LineWidth',2)
xlabel("Supersaturation",'Interpreter','latex','FontSize',10)
ylabel("Primary nucleation rate",'Interpreter','latex','FontSize',10)
xlim([1,1e5])
xticks([1,1e1,1e2,1e3,1e4,1e5])
yticks([1,1e5,1e10,1e15,1e20,1e25])
ylim([1,1e26])
box on
set(gcf,'color','w');
set(gcf,'units','inches','position',[1,1,3.25,2.25])