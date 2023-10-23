function [loggamma12_pure,loggamma14_pure,loggamma32_pure] = pureElectrolyteActivity(catConcs,anConcs) % Na, Cl, OH, Mg

A_gamma=0.511; % valid @ T=25 C

% ions are: Mg Na  indexes are: 1 3
%           OH Cl               2 4

% Electrolytes can be: Mg(OH)2 NaOH
%                      MgCl2   -
z_cations=[2 1];
z_anions=[1 1];

% Related to these indexes: 12 32
%                           14 -

B_=BCalculation;
loggamma_pure=zeros(2);

I_s=ionicStrength(catConcs(1),anConcs(1),anConcs(2),catConcs(2));
sqrt_Is=sqrt(I_s);

% gamma_pure
for j=1:2
    for i=1:2
        B=B_(i,j);
        z_c=z_cations(j); z_a=z_anions(i);
        loggamma_pure(i,j)=-A_gamma*(z_c*z_a)*sqrt_Is/(1+sqrt_Is)+(0.06+0.6*B)*((z_c*z_a))*I_s/(1+1.5*I_s/(z_c*z_a))^2+B*I_s;
    end
end
loggamma12_pure=loggamma_pure(1); loggamma14_pure=loggamma_pure(2); loggamma32_pure=loggamma_pure(3);
end