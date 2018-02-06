% Approximate the real specific acoustical impedance (Z_specific=Z/Z_0) 
% by a uniform absorption degree with the method proposed by
% Mommertz, E. (1996): Untersuchung akustischer Wandeigenschaften und Modellierung
% der Schallrückwürfe in der binauralen Raumsimulation. Dissertation. RWTH
% Aachen, Fakultät für Elektrotechnik und Informationstechnik, S. 122. 

function Z_specific = conversion_alpha_to_Z(alpha_s)

func = @(zeta) (8./(zeta.^2).*(1+zeta-1./(1+zeta)-2.*log(1+zeta))-alpha_s)^2; 
opts = optimset('MaxFunEvals',10000);

% Z_specific_all=0;

for i=0:1000
[Z_specific,error] = fminsearch (func,i);


% break when acoustical impedance is bigger than 1 (bigger than acoustical impedance of air)
if Z_specific>1
    break
end

end

if i==1000;
   msgbox('acoustical impedance cannot be modelled')
end

end

