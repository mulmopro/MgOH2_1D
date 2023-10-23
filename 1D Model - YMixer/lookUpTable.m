function [MgOH2_sol] = lookUpTable(M_MgCl2_,M_NaOH_,variance,s)
% Segregation coefficient matrix construction
N=171;                                           % Nodes number

for j=1:N
    Is=linspace(0,1,N);                          % segregation coefficient
end

% bPDF function A=MgCl2=0, B=NaOH=1;             % convention used
zeta_s = 2*M_MgCl2_/(2*M_MgCl2_ + M_NaOH_);
MgCl2=zeros(1,N);
NaOH=zeros(1,N);
MgOH2=zeros(1,N);

for i=2:N-1
    I_s=Is(i);                                   % segregation coefficient
    ni=s.alpha*(1/I_s-1);
    w=(1-s.alpha)*(1/I_s-1);
    beta=gamma(ni)*gamma(w)/gamma(ni+w);
    beta1=gamma(ni+1)*gamma(w)/gamma(ni+w+1);
    NaOH(i)=(1/(beta-beta*zeta_s)*(gamma(ni+1)*gamma(w)/gamma(ni+w+1)-beta1*betainc(zeta_s,ni+1,w)-zeta_s*beta*betainc(zeta_s,ni,w,"upper")))*M_NaOH_;
    MgCl2(i)=(1/beta*(beta*betainc(zeta_s,ni,w)-1/zeta_s*beta1*betainc(zeta_s,ni+1,w)))*M_MgCl2_;
    MgOH2(i)=(M_NaOH_/beta*beta1*betainc(zeta_s,ni+1,w)+M_NaOH_*zeta_s/(1-zeta_s)/beta*(betainc(zeta_s,ni,w,'upper')*beta-betainc(zeta_s,ni+1,w,'upper')*beta1));
end

for i=1:length(Is)
    if Is(i)==1
        NaOH(i)=M_NaOH_*s.alpha;
        MgCl2(i)=M_MgCl2_*(1-s.alpha);
        MgOH2(i)=0;
    elseif Is(i)==0 && s.alpha<=zeta_s
        NaOH(i)=0;
        MgCl2(i)=(1-s.alpha/zeta_s)*M_MgCl2_;
        MgOH2(i)=M_NaOH_*s.alpha;
    elseif Is(i)==0 && s.alpha>zeta_s
        NaOH(i)=(s.alpha-zeta_s)/(1-zeta_s)*M_NaOH_;
        MgCl2(i)=0;
        MgOH2(i)=M_NaOH_*zeta_s*(1-s.alpha)/(1-zeta_s);
    end
end

MgOH2=spline(Is,MgOH2,linspace(0,1,10000));

OF_Is=variance/(s.alpha*(1-s.alpha));

% Segregation coefficient evaluation
Is=linspace(0,1,length(MgOH2));
MgOH2_sol=interp1(Is,MgOH2,OF_Is);
end