% Drag force
function Fd = theory_drag(ref,parameters)

switch ref
    case 'Stokes-Sublayer-Smooth'  %Soltani(1994)
        rho = parameters.rho;
        Dp = parameters.Dp;
        ustar = parameters.ustar;
        C = parameters.C;
        % Eq. 24
        Fd = 2.9*pi*rho*Dp^2*ustar^2/C;
        return
    case 'Stokes-Burst-Smooth'    %Soltani(1994)
        rho = parameters.rho;
        Dp = parameters.Dp;
        ustar = parameters.ustar;
        C = parameters.C;
        f = 1.7009;
        % Eq. 33 Goldasteh (Neglect nonlinearity)
        Fd = 2.58*pi*f*rho*Dp^2*ustar^2/C;  % (It is correct)
        return
     % double check the following eqs. from soltani 1995 it seems not correct
     case 'Stokes-Sublayer-Rough'
        rho = parameters.rho;
        Dp = parameters.Dp;
        ustar = parameters.ustar;
        C = parameters.C;
        L = parameters.L;
        Fd = 5.8*pi*rho*Dp*ustar^2*L/C;
        return
    case 'Stokes-Burst-Rough'
        rho = parameters.rho;
        nu = parameters.nu;
        Dp = parameters.Dp;
        ustar = parameters.ustar;
        C = parameters.C;
        L = parameters.L;
        Fd = 5.8*pi*rho*Dp*ustar^2*L/C*(1.72 + .1*ustar*L/nu);
        return   
    case 'Nonlinear-Sublayer'   %Goldasteh(2013) ?????????
        % Here gamma = 1.1365 (Burst/inrush intensity factor)
        rho = parameters.rho;
        nu = parameters.nu;
        Dp = parameters.Dp;
        ustar = parameters.ustar;
        C = parameters.C;
        f = 1.7009; % Correction factor for the wall effect
        Res = .56825*(Dp*ustar/nu)^2;   % Particle Reynolds number for the sublayer model (Eq. 28)
        % Eq. 29
        Fd = (1.70475*pi*f*rho*Dp^2*ustar^2/C) * (1 + .15*Res^.678);
        return
    case 'Nonlinear-Burst'   %Goldasteh(2013) ?????????
        % Here gamma = 1.72 (Burst/inrush intensity factor)
        rho = parameters.rho;
        nu = parameters.nu;
        Dp = parameters.Dp;
        ustar = parameters.ustar;
        C = parameters.C;
        f = 1.7009; % Correction factor for the wall effect
        Res = .86*(Dp*ustar/nu)^2;   % Particle Reynolds number for the sublayer model (Eq. 28)
        % Eq. 33
        Fd = (2.58*pi*f*rho*Dp^2*ustar^2/C) * (1 + .15*Res^.678);
        return
end %switch
