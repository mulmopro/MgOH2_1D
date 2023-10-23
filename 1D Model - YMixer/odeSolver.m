function [dydt] = odeSolver(t,y,s)

% disp("Current time: "+num2str(t,'%.2e')+" s")
dmomdt=zeros(2*s.N,1);
dconcdt=zeros(2,1);

if t<=s.tVar(end)
    epsilon_ = epsilon(t,s);
    variance_int=interp1(s.tVar,s.variance,t);
else
    epsilon_=0;
    variance_int=s.variance(end);
end

if s.N==2
    moments=[y(1) y(2) y(3) y(4)];
    % ion order Mg, OH
    concs=[y(5) y(6)];
    % ion order Na, Mg - Cl, OH
    [gammapm] = electrolyteActivity([s.M_NaOH*s.alpha, y(5)],[2*s.M_MgCl2*(1-s.alpha), y(6)],s);
else
    moments=[y(1) y(2) y(3) y(4) y(5) y(6)];
    % ion order Mg, OH
    concs=[y(7) y(8)];
    % ion order Na, Mg - Cl, OH
    [gammapm] = electrolyteActivity([s.M_NaOH*s.alpha, y(7)],[2*s.M_MgCl2*(1-s.alpha), y(8)],s);
end

if y(1)>0 && y(2)>0
    [w,L]=whAlgAdaptive(moments);
else
    w=zeros(1,s.N);
    L=zeros(1,s.N)+eps;
end

MgOH2_sol=lookUpTable(concs(1),concs(2),variance_int,s);

if MgOH2_sol>=power(s.k_sp/(4*gammapm^3),1/3)
    S=(4*gammapm^3*power(MgOH2_sol,3)-s.k_sp)/s.k_sp;
else
    S=0;
end

J=s.A_*exp(-s.B_/power(log(S+1),2))+s.A_het*exp(-s.B_het/power(log(S+1),2));
if s.diffContr==false
    G=s.kg_*S^s.g_;
else
    for i=1:s.N
        G(i)=diffusionControlled(t,L(i),gammapm,S);
    end
end
% aggregation evaluation

betaKern=zeros(s.N,s.N);
for i=1:s.N
    for j=1:s.N
        if t<=s.tVar(end) && S>0
        betaCol= s.c_hd*sqrt(epsilon_/s.nu)*power((L(i)+L(j)),3) + 2*s.kb*s.T/(3*s.mu)*(L(i)/L(j)+L(j)/L(i)+2); %  
        betaKern(i,j)=aggregEff(L(i),L(j),epsilon_,mean(G),S,s)*betaCol*10^(s.C_adjust);
        else
        betaKern(i,j)=0;
        end
    end
end

aggregRate=zeros(2*s.N,1);
for k=1:2*s.N
    aggregRateSum=0;
    for i=1:s.N
        for j=1:s.N
            aggregRateSum=aggregRateSum+0.5*w(i)*(w(j)*power(L(i)^3+L(j)^3,(k-1)/3))*betaKern(i,j)-w(i)*L(i)^(k-1)*(betaKern(i,j)*w(j));
        end
    end
    aggregRate(k)=aggregRateSum;
end

if s.N==2
    if s.aggregation==false
        dmomdt(1)=J;
        dmomdt(2)=J*(s.nucleateSize)+sum(G.*w);
        dmomdt(3)=J*(s.nucleateSize)^2+2*sum(G.*w.*(L));
        dmomdt(4)=J*(s.nucleateSize)^3+3*sum(G.*w.*(L.^2));
    else
        dmomdt(1)=J+aggregRate(1);
        dmomdt(2)=J*(s.nucleateSize)+sum(G.*w)+aggregRate(2);
        dmomdt(3)=J*(s.nucleateSize)^2+2*sum(G.*w.*(L))+aggregRate(3);
        dmomdt(4)=J*(s.nucleateSize)^3+3*sum(G.*w.*(L.^2))+aggregRate(4);
    end
else
    if s.aggregation==false
        dmomdt(1)=J;
        dmomdt(2)=J*(s.nucleateSize)+G*y(1);
        dmomdt(3)=J*(s.nucleateSize)^2+2*G*y(2);
        dmomdt(4)=J*(s.nucleateSize)^3+3*G*y(3);
        dmomdt(5)=J*(s.nucleateSize)^4+4*G*y(4);
        dmomdt(6)=J*(s.nucleateSize)^5+5*G*y(5);
    else
        dmomdt(1)=J+aggregRate(1);
        dmomdt(2)=J*(s.nucleateSize)+sum(G.*w)+aggregRate(2);
        dmomdt(3)=J*(s.nucleateSize)^2+2*sum(G.*w.*(L))+aggregRate(3);
        dmomdt(4)=J*(s.nucleateSize)^3+3*sum(G.*w.*(L.^2))+aggregRate(4);
        dmomdt(5)=J*(s.nucleateSize)^4+4*sum(G.*w.*(L.^3))+aggregRate(5);
        dmomdt(6)=J*(s.nucleateSize)^5+5*sum(G.*w.*(L.^4))+aggregRate(6);
    end
end

dconcdt(1)=-s.kv*dmomdt(4)*s.rhoCry/s.PMCry;
dconcdt(2)=2*dconcdt(1);

dydt=cat(1,dmomdt,dconcdt);
end