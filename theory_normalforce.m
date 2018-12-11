function Fn = theory_normalforce(ref,parameters)
switch ref
    case 'With lift/gravitational forces'
        Dp = parameters.Dp;
        Fpo = parameters.Fpo;
        rho_p = parameters.rho_p;
        Fmg = rho_p * 4/3*pi*(Dp/2)^3 * 9.81;
        Fl = parameters.Fl;
    case 'No lift/gravitational forces'
        Fpo = parameters.Fpo;
        Fmg = 0;
        Fl = 0;
    case 'No gravitational force'
        Fpo = parameters.Fpo;
        Fmg = 0;
        Fl = parameters.Fl;
end %switch
Fn = Fpo + Fmg - Fl;