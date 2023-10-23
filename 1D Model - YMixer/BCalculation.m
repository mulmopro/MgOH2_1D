function [B] = BCalculation()
% ions are: Mg Na  indexes are: 1 3
%           OH Cl               2 4
% Electrolytes can be: Mg(OH)2 NaOH
%                      MgCl2   -


% B_Mg, delta_Mg
B_Mg=0.057; delta_Mg=0.157;                        % kg/mol, sqrt(kg/mol)
% B_OH, delta_OH
B_OH=0.076; delta_OH=-1;                           % kg/mol, sqrt(kg/mol)

B=[B_Mg+B_OH+delta_Mg*delta_OH, 0.0747;0.1129,0];

end