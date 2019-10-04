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
        
    case 'Rumpf-Rabinovich'
        % Rump-Rabinovich adhesion force (Modified Rumpf model)
        % Referece: Rabinovich et al. (2000) Eq. (7) [first paper]
        % here H0=z0, R=Dp/2, sigma: rms
        R = Dp/2;
        H0 = parameters.z0;
        sigma = parameters.sigma;
        if isnan(parameters.A)
            A = Wa * (12*pi*H0^2);
        else
            A = parameters.A;
        end %if
        Fpo = (A*R/(6*H0^2)) .* ((1+R./(1.48.*sigma)).^(-1) + (1+1.48.*sigma/H0).^(-2));
        return
        
    case 'Rabinovich'
        R = Dp/2;
        H0 = parameters.z0;
        sigma = parameters.sigma;
        if isnan(parameters.A)
            A = Wa * (12*pi*H0^2);
        else
            A = parameters.A;
        end %if
        lambda = parameters.lambda;
        
        k1 = 1.817;
        Fpo = (A*R/(6*H0^2)) .* ( (1+32*k1*R.*sigma./lambda.^2).^(-1) + (1+k1*sigma/H0).^(-2) );
        return
        
end %switch
