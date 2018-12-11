% Pull-off force
% Wa : Thermodynamic work of adhesion (Surface energy per unit area)
function Fpo = theory_pulloff(model,parameters)
Dp = parameters.Dp;
Wa = parameters.Wa;
switch model
    case 'JKR'
        Fpo = 3/4*pi*Wa*Dp;
    case 'DMT'
        Fpo = pi*Wa*Dp;
    case 'TPL'
        A = parameters.A;
        z0 = parameters.z0;
        K = parameters.K;
        V = (25*A^2*Dp/(288*z0^7*K^2))^(1/3);
        Fo = pi*Wa*Dp;
        Fpo = Fo*(.5*exp(.124*(V-.01)^.439) + .2*V);
    case 'JKR-asperity'
        beta = parameters.betap * Dp;
        Fpo = 3/2*pi*Wa*beta;
    case 'JKR-Rough'
        a = parameters.a;
        N = parameters.N;
        fpo = parameters.fpo;
        Deltac = parameters.Deltac;
        Fpo = pi*a^2*N*fpo*exp(-.6/Deltac^2);
    case 'JKR-Rough-Iman'
        Deltac = parameters.Deltac;
        Fpo = .0029/4*Dp*Wa * (1.5*pi^2*exp(-.6/Deltac^2))^3;
    case 'JKR-bump'
        K = parameters.K;
%         N = parameters.N;
        sigma = parameters.sigma;
%         beta = parameters.betap * Dp;
        Deltac = parameters.Deltac;
        n_u = parameters.n_u;
        n_b = parameters.n_b;
        Nbump = parameters.Nbump;
        Fpo = (Dp/(n_u*n_b*sqrt(Nbump)*K))^2 * (3/2*(0.1/sigma)*pi^2*Wa*exp(-.6/Deltac^2))^3;
    case 'JKR-bump-new' % It is not validate yet!!!!! (Iman's Thesis)
        Deltac = parameters.Deltac;
        n_u = parameters.n_u;
        n_b = parameters.n_b;
        Nbump = parameters.Nbump;
        Fpo = (Dp/(n_u*n_b*sqrt(Nbump)))*.00145*(1.5*pi^2*Wa^3*exp(-.6/Deltac^2))^3;
end %switch
