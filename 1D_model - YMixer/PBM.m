function [tSol,ySol] = PBM(s)

% ODESolver
if s.N==2
    options=odeset;
    [tSol,ySol]=ode15s(@odeSolver,[s.tVar(1),s.tend],[0 0 0 0 s.M_MgCl2*(1-s.alpha) s.M_NaOH*s.alpha],options,s);
else
    options=odeset('AbsTol',1e-16);
    [tSol,ySol]=ode15s(@odeSolver,[s.tVar(1),s.tend],[0 0 0 0 0 0 s.M_MgCl2*(1-s.alpha) s.M_NaOH*s.alpha],options,s);
end

end
