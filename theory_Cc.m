% Cunningham correction factor
% Dp     : Particle Diameter (m)
% lambda : Gas mean free path (m)
function Cc = theory_Cc(Dp, lambda,const)
if (isempty(const) == 1)
%   Constants from Prof. Ahmadi's notes
%   Refer to Fuchs[39] and Friedlander[40] which are References of Soltani (1994)
    alpha = 1.257;
    beta = 0.4;
    gamma = 1.1;
else
    alpha = const.alpha;
    beta = const.beta;
    gamma = const.gamma;
end %if
Kn = 2*lambda./Dp;   % Knudsen Number
Cc = 1 + Kn.*(alpha + beta.*exp(-gamma./Kn));
