% The moment of surface stresses
function M = theory_momentp(ref,parameters)
switch ref
    case 'Sublayer-Smooth'
        rho = parameters.rho;
        ustar = parameters.ustar;
        Dp = parameters.Dp;
        C = parameters.C;
        fm = .943993;
        % Eq. 30 (Soltani)
        % Eq. 31 (Goldasteh)
        M = 1.1365*fm*pi*rho*ustar^2*Dp^3/C;
        return
    case 'Burst-Smooth'
        rho = parameters.rho;
        ustar = parameters.ustar;
        Dp = parameters.Dp;
        C = parameters.C;
        fm = .943993;
        % Eq. 35 (Goldasteh)
        M = 1.72*pi*rho*fm*Dp^3*ustar^2/C; % (It is correct)
        return
    % Double check the following Eqs. for rough (from soltani 1995 it seems
    % not correct)
    case 'Sublayer-Rough'
        rho = parameters.rho;
        ustar = parameters.ustar;
        Dp = parameters.Dp;
        C = parameters.C;
        L = parameters.L;
%         M = 2.14*pi*rho*ustar^2*Dp^2*L/C;  %Soltani(1995)
        fm = .943993;
        gamma = 1.1365; %sublayer
        M = gamma*pi*rho*fm*Dp^3*ustar^2/C;
        return
    case 'Burst-Rough'
        rho = parameters.rho;
        nu = parameters.nu;
        ustar = parameters.ustar;
        Dp = parameters.Dp;
        C = parameters.C;
        L = parameters.L;
%         M = 2.14*pi*rho*Dp^2*ustar^2*L/C*(1.72 + .1*ustar*L/nu); % Soltani(1995)
        fm = .943993;
        gamma = 1.72;
        M = gamma*pi*rho*fm*Dp^3*ustar^2/C;
        return
end %switch
