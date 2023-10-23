function [gamma] = electrolyteActivity(catConcs,anConcs,s)  % Na, Cl, OH, Mg

if s.activity==true
    A_gamma=0.511;
    rho_solvent=997;                                        %kg/m^3

    y_ac=zeros(1,2);
    x_ca=zeros(1,2);
    z_meanSquare=[9 9; 9 4]; % first row is for Y (21-41), second row is for X (12-32)

    [loggamma12_pure,loggamma14_pure,loggamma32_pure] = pureElectrolyteActivity(catConcs,anConcs);
    I_s=ionicStrength(catConcs(1),anConcs(1),anConcs(2),catConcs(2));
    sqrt_Is=sqrt(I_s);

    y_ac(1)=z_meanSquare(1,1)*(anConcs(2)*1000/rho_solvent)/(4*I_s);
    y_ac(2)=z_meanSquare(1,2)*(anConcs(1)*1000/rho_solvent)/(4*I_s);

    x_ca(1)=z_meanSquare(2,1)*(catConcs(2)*1000/rho_solvent)/(4*I_s);
    x_ca(2)=z_meanSquare(2,2)*(catConcs(1)*1000/rho_solvent)/(4*I_s);

    % "Mixing rules" implementation
    f_c=y_ac(1)*loggamma12_pure+y_ac(2)*loggamma14_pure+A_gamma*sqrt_Is/(1+sqrt_Is)*(2*y_ac(1)+2*y_ac(2));
    f_a=x_ca(1)*loggamma12_pure+x_ca(2)*loggamma32_pure+A_gamma*sqrt_Is/(1+sqrt_Is)*(2*x_ca(1)+x_ca(2));

    % it's for Mg(OH)2
    z_c=2;
    z_a=1;

    nu_c=z_a;
    nu_a=z_c;
    nu=nu_c+nu_a;

    loggamma=-A_gamma*(z_c*z_a)*sqrt_Is/(1+sqrt_Is)+nu_c*f_c/nu+nu_a*f_a/nu;
    gamma=10^loggamma;
else
    gamma=1;
end

end