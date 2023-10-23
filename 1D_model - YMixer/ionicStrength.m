function [Im] = ionicStrength(Na_,Cl_,OH_tot_,Mg_tot_)
    
        rho=997;                                                                               % g/l
        z_Na=1; z_Cl=-1; z_Mg=2; z_OH=-1;
        m_Na=Na_*1000/rho; m_Cl=Cl_*1000/rho; m_Mg=Mg_tot_*1000/rho; m_OH=OH_tot_*1000/rho;    % mol/kg_solv
        
        % Ionic strength calculation
        Im=0.5*(m_Na*power(z_Na,2)+m_Cl*power(z_Cl,2)+m_Mg*power(z_Mg,2)+m_OH*power(z_OH,2));  % mol/kg_solv 
end

