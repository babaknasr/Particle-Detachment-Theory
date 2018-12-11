% Resistance Moment
function rmoment = theory_rmoment(ref,parameters)
% if ~isfield(parameters,'assumption')
%     parameters.assumption.gravityforce = 0;
% else
%     if ~isfield(parameters.assumption,'gravityforce')
%         parameters.assumption.gravityforce = 0;
%     end %if
% end %if
if (parameters.assumption.gravityforce == 1)
    Dp = parameters.Dp;
    rho_p = parameters.rho_p;
    Fmg = rho_p * 4/3*pi*(Dp/2)^3 * 9.81;
else
    Fmg = 0;
end %if
switch ref
    case 'Detached moment-Smooth'
%         Dp = parameters.Dp;
        a = parameters.a;
        Fpo = parameters.Fpo;
%         rho_p = parameters.rho_p;
%         Fmg = rho_p * 4/3*pi*(Dp/2)^3 * 9.81;
        rmoment = (Fpo + Fmg)*a;
        return
    case 'Detached moment-No gravitational force-Smooth'
        a = parameters.a;
        Fpo = parameters.Fpo;
        Fmg = 0;
        rmoment = (Fpo + Fmg)*a;
        return
    case 'MmaxJKR-Smooth'
        Dp = parameters.Dp;
        Wa = parameters.Wa;
        K = parameters.K;
        a = parameters.a;
%         rho_p = parameters.rho_p;
        Mmax = 2.70716 * Wa^(4/3)*Dp^(5/3)/K^(1/3);
%         Fmg = rho_p * 4/3*pi*(Dp/2)^3 * 9.81;
        rmoment = Mmax + Fmg*a;
        return
    case 'MmaxDMT-Smooth'
        Dp = parameters.Dp;
        Wa = parameters.Wa;
        K = parameters.K;
        a = parameters.a;
%         rho_p = parameters.rho_p;
        Mmax = 1.7254 * Wa^(4/3)*Dp^(5/3)/K^(1/3);
%         Fmg = rho_p * 4/3*pi*(Dp/2)^3 * 9.81;
        rmoment = Mmax + Fmg*a;
        return
    case 'Detached moment-Rough'
%         Dp = parameters.Dp;
        a = parameters.a;
%         rho_p = parameters.rho_p;
        FM = parameters.FM;
%         Fmg = rho_p * 4/3*pi*(Dp/2)^3 * 9.81;
        rmoment = (FM + Fmg)*a;
        return
end %switch

