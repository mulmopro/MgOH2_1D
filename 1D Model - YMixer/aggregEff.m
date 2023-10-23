
function [Pa] = aggregEff(L1,L2,eps_,G,superSaturation,s)

if s.aggregEfficiency==true
    if superSaturation > 0 && G > 0 && L1 >= 1e-9 && L2 >= 1e-9
        ti=sqrt(s.nu/eps_);                                   % s
        Leq=(L1*L2)/sqrt(L1^2+L2^2-L1*L2);
        Db=Leq*sqrt(s.rhoCry/power(10,s.Ap))*(eps_*s.nu)^0.25;
        
        if L1>=L2
            r_L=L1/L2;
        else
            r_L=L2/L1;
        end
        
        sqrt_r_L=sqrt(r_L^2-1);
        f_lambda=4 * (1 + r_L - sqrt_r_L)/ ((1/3 + r_L - sqrt_r_L)-((r_L - sqrt_r_L)^2)*(2*r_L/3 + sqrt_r_L/3));
        tc=Db/f_lambda/G;
        Pa=exp(-tc/ti);
    else
        Pa=0;
    end
else
        Pa=1;
end

end