function Fmg = theory_gravity(ref,parameters)
Dp = parameters.Dp;
switch ref
    case 'No gravity'
        Fmg = 0;
    case 'With gravity'
        rho_p = prameters.rho_p;
        Fmg = rho_p * 4/3*pi*(Dp/2)^3 * 9.81;
end %switch