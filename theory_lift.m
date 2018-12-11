% Lift force
function Fl = theory_lift(ref,parameters)
% Soltani(1994)
switch ref
    case 'Saffman-Sublayer-Smooth'
        rho = parameters.rho;
        nu = parameters.nu;
        Dp = parameters.Dp;
        ustar = parameters.ustar;
        % Eq. 28
        Fl = .975*rho*nu^2*(Dp*ustar/nu)^3;
        return
    case 'Saffman-Burst-Smooth'
        rho = parameters.rho;
        nu = parameters.nu;
        Dp = parameters.Dp;
        ustar = parameters.ustar;
        Fl = 1.8159*rho*nu^2*(Dp*ustar/nu)^3;
        return
    case 'Saffman-Sublayer-Rough'
        rho = parameters.rho;
        nu = parameters.nu;
        Dp = parameters.Dp;
        ustar = parameters.ustar;
        L = parameters.L;
        Fl = 1.95*rho*Dp^2*ustar^3*L/nu;
        return
    case 'Saffman-Burst-Rough'
        rho = parameters.rho;
        nu = parameters.nu;
        Dp = parameters.Dp;
        ustar = parameters.ustar;
        L = parameters.L;
        Fl = 1.95*rho*Dp^2*ustar^3*L/nu*(1.72 + .1*ustar*L/nu)*(1.72 + .2*ustar*L/nu)^.5;
        return
end %switch


% % from Prof. Ahmadi's notes (Lift force Eq.16 & 17)
% dpplus = Dp_plus(Dp,ustar,nu);
% Flplus = 0;
% if (dpplus > 1.5)
%     Flplus = 4.21*dpplus^2.31;  % Hall(1988)
% end %if
% if ((dpplus > .15) && (dpplus < 1))
%     Flplus = 15.57*dpplus^1.87; % Nieuwstadt (1996)
% end %if
% Fl = Flplus*rho*nu^2;
