function [d]=MgOH2Plot(tSol,ySol,s)

tempo=tSol;
if s.N==2
    m0=ySol(:,1);
    m1=ySol(:,2);
    m2=ySol(:,3);
    m3=ySol(:,4);
    Mgsolution=ySol(:,5);
    OHsolution=ySol(:,6);
else
    m0=ySol(:,1);
    m1=ySol(:,2);
    m2=ySol(:,3);
    m3=ySol(:,4);
    m4=ySol(:,5);
    m5=ySol(:,6);
    Mgsolution=ySol(:,7);
    OHsolution=ySol(:,8);
end

% memory pre-allocation
MgOH2_sol=zeros(length(tempo),1); ionicStrength_=zeros(length(tempo),1); superSaturation=zeros(length(tempo),1); gamma=zeros(length(tempo),1); J=zeros(length(tempo),1); J_het=zeros(length(tempo),1); G=zeros(length(tempo),1); Im=zeros(length(tempo),1); w=zeros(length(tempo),s.N); L=zeros(length(tempo),s.N);

for i=1:length(tempo)
    variance_int=interp1(s.tVar,s.variance,tempo(i));
    if variance_int<0
        variance_int=0;
    end

    MgOH2_sol(i)=lookUpTable(Mgsolution(i),OHsolution(i),variance_int,s);
    ionicStrength_(i)=ionicStrength(s.M_NaOH*s.alpha,s.M_MgCl2*s.alpha,OHsolution(i),Mgsolution(i));
    gammapm = electrolyteActivity([s.M_NaOH*s.alpha, Mgsolution(i)],[s.M_MgCl2*s.alpha, OHsolution(i)],s);
    gamma(i)=gammapm;

    if MgOH2_sol(i)>power(s.k_sp/(4*gamma(i)^3),1/3)
        superSaturation(i)=(gamma(i)^3*4*MgOH2_sol(i)^3-s.k_sp)/s.k_sp;
    else
        superSaturation(i)=0;
    end
    J(i)=s.A_*exp(-s.B_/power(log(superSaturation(i)+1),2))+s.A_het*exp(-s.B_het/power(log(superSaturation(i)+1),2));
    G(i)=s.kg_*superSaturation(i)^s.g_;

    if m0(i)>0 && m1(i)>0
        if s.N==2
            [w(i,:),L(i,:)]=whAlgAdaptive([m0(i) m1(i) m2(i) m3(i)]);
        else
            [w(i,:),L(i,:)]=whAlgAdaptive([m0(i) m1(i) m2(i) m3(i) m4(i) m5(i)]);
        end
    else
        w(i,:)=zeros(1,s.N);
        L(i,:)=zeros(1,s.N)+eps;
    end
end

% Particle size distribution reconstruction
% [particleSizes,logNormDist] = reconstructNSD(m0(end),m1(end),m2(end));

% Results
% if s.N==2
%     m3_teor=(s.PMCry*(s.M_MgCl2*s.alpha-ySol(end,5))/(s.kv*s.rhoCry));
% else
%     m3_teor=(s.PMCry*(s.M_MgCl2*s.alpha-ySol(end,7))/(s.kv*s.rhoCry));
% end
%
% disp("Mass balance error is: "+num2str(abs(m3_teor-m3(end))/m3_teor*100) + " %")
% disp("Final Mg conversion is: " + num2str(((s.M_MgCl2*s.alpha-Mgsolution(end))/s.M_MgCl2/s.alpha)*100)+ " %")
% disp("Final d10 value is: " + num2str(m1(end)/m0(end)*10^6)+ " micrometers")
% disp("Final d21 value is: " + num2str(m2(end)/m1(end)*10^6)+ " micrometers")
% disp("Final d32 value is: " + num2str(m3(end)/m2(end)*10^6)+ " micrometers")
% disp("Final d43 value is: " + num2str(m4(end)/m3(end)*10^6)+ " micrometers")
%
% % Plotting fields
% figure('DefaultAxesFontSize',16)
% subplot(2,s.N,1), plot(tempo,m0,'LineWidth',2.5)
% ylabel('$m_{0},$\ $\frac{number\ of\ particles}{m^{3}_{suspension}}$','Interpreter','latex','FontSize',24)
% subplot(2,s.N,2), plot(tempo,m1,'LineWidth',2.5)
% ylabel('$m_{1},$\ $\frac{m}{m^{3}_{suspension}}$','Interpreter','latex','FontSize',24)
% subplot(2,s.N,3), plot(tempo,m2,'LineWidth',2.5)
% ylabel('$m_{2},$\ $\frac{m^{2}}{m^{3}_{suspension}}$','Interpreter','latex','FontSize',24)
% if s.N==2
%     xlabel('time, s','FontSize',16)
% end
% subplot(2,s.N,4), plot(tempo,m3,'LineWidth',2.5)
% ylabel('$m_{3},$\ $\frac{m^{3}}{m^{3}_{suspension}}$','Interpreter','latex','FontSize',24)
% xlabel('time, s','FontSize',16)
% if s.N==3
%     subplot(2,s.N,5), plot(tempo,m4,'LineWidth',2)
%     ylabel('$m_{4},$\ $\frac{m^{4}}{m^{3}_{suspension}}$','Interpreter','latex','FontSize',24)
%     xlabel('time, s','FontSize',16)
%     subplot(2,s.N,6), plot(tempo,m5,'LineWidth',2)
%     ylabel('$m_{5},$\ $\frac{m^{5}}{m^{3}_{suspension}}$','Interpreter','latex','FontSize',24)
%     xlabel('time, s','FontSize',16)
% end
% set(gcf,'color','w');
%
% figure('DefaultAxesFontSize',10)
% semilogy(1./((log(superSaturation+1)).^2),J,'-k','LineWidth',2), hold on, semilogy(1./((log(superSaturation+1)).^2),J2,'--k','LineWidth',2)
% xlabel("$\frac{\rm{1}}{\rm{ln^{2}(S+1)}}$","Interpreter","latex",'FontSize',10)
% ylabel({'nucleation rate'; '$\frac{\rm{number\ of\ particles}}{\rm{m^{3}s}}$'},'Interpreter','latex','FontSize',10)
% xlim([0 0.17])
% xticks([1e-2,0.05,0.1,0.15])
% box on;
% set(gcf,'color','w');
% set(gcf,'units','inches','position',[1,1,3.25,3.25])
% legend("homogeneous","heterogeneous","FontSize",12,'Location','best')
%
% subplot(2,1,2), semilogx(tempo,G,'LineWidth',2)
% ylabel('growth, $\frac{m}{s}$','Interpreter','latex','FontSize',20)
% xlabel('time, s','FontSize',20)
%
% set(gcf,'color','w');
%
% figure('DefaultAxesFontSize',16)
% plot(tempo,superSaturation,'LineWidth',2)
% ylabel('superSaturation','FontSize',20)
% xlabel('time, s','FontSize',20)
% set(gcf,'color','w');
%
% figure('DefaultAxesFontSize',16)
% subplot(2,2,1), plot(tempo,m1./m0*10^9,'LineWidth',2)
% ylabel('d_{10}, nm','FontSize',25)
% subplot(2,2,2), plot(tempo,m2./m1*10^9,'LineWidth',2)
% ylabel('d_{21}, nm','FontSize',25)
% subplot(2,2,3), plot(tempo,m3./m2*10^9,'LineWidth',2)
% ylabel('d_{32}, nm','FontSize',25)
% xlabel('time, s','FontSize',25)
%
% if s.N==3
%     subplot(2,2,4), plot(tempo,m4./m3*10^9,'LineWidth',2)
%     ylabel('d_{43}, nm','FontSize',25)
%     xlabel('time, s','FontSize',25)
% end
%
% set(gcf,'color','w');
%
% figure('DefaultAxesFontSize',16)
% subplot(1,2,1), plot(tempo,Mgsolution,'LineWidth',2)
% legend('(c_{Mg^{+2}})')
% xlabel('time, s','FontSize',18)
% ylabel('Concentration, M','Interpreter','latex','FontSize',18)
% subplot(1,2,2), plot(tempo,OHsolution,'LineWidth',2)
% legend('(c_{OH^{-}})')
% xlabel('time, s','FontSize',18)
% ylabel('Concentration, M','FontSize',18)
%
% figure('DefaultAxesFontSize',16)
% subplot(1,s.N,1)
% yyaxis right
% plot(tempo,w(:,1),'LineWidth',2)
% ylabel('w_{0}')
% yyaxis left
% plot(tempo,L(:,1),'LineWidth',2)
% ylabel('L_{0}')
% xlabel('time, s','FontSize',18)
%
% subplot(1,s.N,2)
% yyaxis right
% plot(tempo,w(:,2),'LineWidth',2)
% ylabel('w_{1}')
% yyaxis left
% plot(tempo,L(:,2),'LineWidth',2)
% ylabel('L_{1}')
% xlabel('time, s','FontSize',18)
%
% if s.N==3
%     subplot(1,s.N,3)
%     yyaxis right
%     plot(tempo,w(:,3),'LineWidth',2)
%     ylabel('w_{2}')
%     yyaxis left
%     plot(tempo,L(:,3),'LineWidth',2)
%     ylabel('L_{2}')
%     xlabel('time, s','FontSize',18)
% end
% set(gcf,'color','w');
%
% figure('DefaultAxesFontSize',16)
% yyaxis left
% plot(tempo,gamma,'LineWidth',2)
% xlabel('time, s')
% ylabel('$\gamma_{\pm}$','Interpreter','latex','FontSize',18)
% yyaxis right
% plot(tempo,ionicStrength_,'LineWidth',2)
% ylabel('Ionic strength, $\frac{mol}{kg_{solvent}}$','Interpreter','latex','FontSize',18)
% set(gcf,'color','w');
%
% figure('DefaultAxesFontSize',16)
% plot(tempo,Nasolution,'LineWidth',2), hold on, plot(tempo,Clsolution,'LineWidth',2)
% xlabel('time, s')
% ylabel('$Na^{+},\ M$','Interpreter','latex','FontSize',18)
%
% figure('DefaultAxesFontSize',16)
% plot(tempo,V,'LineWidth',2)
% xlabel('time, s')
% ylabel('$Volume,\ l$','Interpreter','latex','FontSize',18)
%
% figure('DefaultAxesFontSize',16)
% % semilogx(particleSizes(1:end-1),n,'LineWidth',2), hold on,
% semilogx(particleSizes,logNormDist,'LineWidth',2),
% fprintf('\n')
% disp("m0 is: "+num2str(sum(n)))
% xlabel('particle size, \mum','FontSize',18)
% ylabel('q0, 1/\mum','FontSize',18)
% legend('Simulation')
% xlim([0.08,20])
% set(gcf,'color','w');
d = [m1(end)/m0(end), m2(end)/m1(end), m3(end)/m2(end), m4(end)/m3(end)];
end