function G = diffusionControlled(t,L,gammapm,S)
global nu_ T_ k_sp_ PMCry_ rhoCry_ nucleateSize tVar

if S>0 && L>nucleateSize
    D=1e-9;                            % m^2/s
    interTension=0.095;                % J/m^2
    R=8.314;                           % J/(mol*K)
    C_sInf=(k_sp_/(4*gammapm))^(1/3);
    v_mol=PMCry_/rhoCry_*0.001;        %m^3/mol
    expFactor=2*interTension*v_mol./(R*T_*L);
    C_sCurv=C_sInf*exp(expFactor);
    if t<=tVar(end)
        Re=epsilon(t).^(1/3)*L.^(4/3)./nu_;
    else
        Re=0;
    end
    Sc=nu_./D;
    Sh=2+0.52*Re.^0.52*Sc.^(1/3);
    
    G=2*(Sh*D./L).*v_mol.*C_sCurv*S;
else
    G=0;
end
end

