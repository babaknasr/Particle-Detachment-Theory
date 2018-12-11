% Contact radius at the moment of separation
function a = theory_cradius(model,parameters)
Dp = parameters.Dp;
Wa = parameters.Wa;
K = parameters.K;
switch model
    case 'JKR-a0'
        a = (3*pi*Wa*Dp^2/(2*K))^(1/3);
    case 'JKR'
        a = (3*pi*Wa*Dp^2/(8*K))^(1/3);
    case 'JKR-Mmax'
        a = .9205 * (3*pi*Wa*Dp^2/(4*K))^(1/3);
    
    case 'DMT'
        a = 0;
    case 'DMT-Mmax'
        a = .5868 * (3*pi*Wa*Dp^2/(4*K))^(1/3);
    case 'TPL'
        A = parameters.A;
        z0 = parameters.z0;
        V = (25*A^2*Dp/(288*z0^7*K^2))^(1/3);
        if (V <= 1.6)
            K20 = .885*(exp(.8*V^.5) - 1);
        else
            K20 = .735*V^.178 + .52*V;
        end %if
        a = Dp*sqrt(K20*z0/(2*Dp));
    case 'JKR-Rough'
        N = parameters.N;
        fpo = parameters.fpo;
        Deltac = parameters.Deltac;
        a = pi*N*fpo*Dp*exp(-.6/Deltac^2)/(2*K);
    case 'JKR-Rough-Iman'
        Deltac = parameters.Deltac;
        a = 3/4*pi^2 * (.0029*Wa*Dp^2/K)^(1/3) * exp(-.6/Deltac^2);
    case 'JKR-Rough-bump'
        N = parameters.N;
        fpo = parameters.fpo;
        bumpr = parameters.bumpr;
        beta = parameters.betap * Dp;
        Deltac = parameters.Deltac;
        a = pi*N*fpo*bumpr*exp(-.6/Deltac^2)/K;
end %switch
