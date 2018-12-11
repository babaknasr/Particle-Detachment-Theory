function ustar = theory_vc_soltani(model,mode,parameters)
global glob_sigma glob_N glob_gamma

ustar = -1;
rho = parameters.rho;
nu = parameters.nu;
z0 = parameters.z0;
Wa = parameters.Wa;
A = parameters.A;
K = parameters.K;
Dp = parameters.Dp;
C = parameters.C;
mus = parameters.mus;
rho_p = parameters.rho_p;
us00 = 0; %!!!

switch model
% *********** Sublayer ***********
% *********** Sublayer  - Rolling     ***************
    case 'TPL-Rolling-Sublayer-Smooth'    
        V = (25*A^2*Dp/(288*z0^7*K^2))^(1/3);
        if (V <= 1.6)
            K20 = .885*(exp(.8*V^.5) - 1);
        else
            K20 = .735*V^.178 + .52*V;
        end %if
        switch mode
            case 'No gravitational force'
                us0 = abs(((Wa*pi*(.5*exp(.124*(V-.01)^.439) + .2*V) - .975*rho/nu*Dp^2*us00) * .28*C/(rho*pi)*sqrt(K20*z0/Dp^3))^.5);
                fu = @(us) us - ((Wa*pi*(.5*exp(.124*(V-.01)^.439) + .2*V) - .975*rho/nu*Dp^2*us^3) * .28*C/(rho*pi)*sqrt(K20*z0/Dp^3))^.5;
                options = optimset('TolFun',1e-10,'Tolx',1e-10);
                ustar = fzero(fu,us0,options);
            case 'No lift/gravitational forces'
                ustar = (.28*C*Wa/rho*(.5*exp(.124*(V-.01)^.439) + .2*V)*sqrt(K20*z0/Dp^3))^(1/2);
            case 'No hyd. torque'
                ustar = (.55*C*Wa/rho*(.5*exp(.124*(V-.01)^.439) + .2*V)*sqrt(K20*z0/Dp^3))^(1/2);
        end %switch mode
    case 'JKR-Rolling-Sublayer-Smooth'
        switch mode
            case 'No lift/gravitational forces'
                ustar = .46*(C^(3/2)*Wa^2*pi^.5/(rho^(3/2)*Dp^2*K^.5))^(1/3);
        end %switch mode
%*********** Sublayer - Sliding **********        
    case 'TPL-Sliding-Sublayer-Smooth'
        V = (25*A^2*Dp/(288*z0^7*K^2))^(1/3);
        switch mode
            case 'No lift/gravitational forces' %??????
                ustar = (.34*C*mus*Wa/(rho*Dp)*(.5*exp(.124*(V-.01)^.439) + .2*V))^.5;
        end %switch
    case 'JKR-Sliding-Sublayer-Smooth'
        switch mode
            case 'No lift/gravitational forces' %??????
                ustar = .5*(C*mus*Wa/(rho*Dp))^.5;
        end %switch
%************* Burst/Inrush
%************* Burst/Inrush - Rolling ***********
    case 'TPL-Rolling-Burst-Smooth'
        V = (25*A^2*Dp/(288*z0^7*K^2))^(1/3);
        if (V <= 1.6)
            K20 = .885*(exp(.8*V^.5) - 1);
        else
            K20 = .735*V^.178 + .52*V;
        end %if
        switch mode
            case 'No gravitational force' %?????
                us0 = ((.5*exp(.124*(V-.01)^.439) + .2*V) * .32*C*Wa/rho * sqrt(K20*z0/Dp^3)/(2.43 + .07*Dp*us00/nu))^.5;
                fu = @(us) us - ((.5*exp(.124*(V-.01)^.439) + .2*V) * .32*C*Wa/rho * sqrt(K20*z0/Dp^3)/(2.43 + .07*Dp*us/nu))^.5;
                options = optimset('TolFun',1e-10,'Tolx',1e-10);
                ustar = fzero(fu,us0,options);
        end %switch
    case 'JKR-Rolling-Burst-Smooth'
        switch mode
            case 'No gravitational force' %??????
                us0 = .5*(C^(3/2)*Wa^2*pi^(1/2)/(rho^(3/2)*K^.5*Dp^2))^(1/3) / (2.43 + .07*Dp*us00/nu)^.5;
                fu = @(us) us - .5*(C^(3/2)*Wa^2*pi^(1/2)/(rho^(3/2)*K^.5*Dp^2))^(1/3) / (2.43 + .07*Dp*us/nu)^.5;
                options = optimset('TolFun',1e-15,'Tolx',1e-15);
                ustar = fzero(fu,us0,options);
            case 'Max moment' % Ahmadi(2007) Stoke's drag
                Mmax = 2.707*(Wa^4*Dp^5/K)^(1/3);
                ustar = sqrt(Mmax/(3.81*pi*rho*Dp^3/C));
                
            case 'Detached moment (Babak)' % Detached moment and nonlinear drag
                gamma = 1.72; %burst
                fD = 1.7009;
                fM = .943993;
                Mad = 3*pi*(3*pi)^(1/3)/8*(Wa^4*Dp^5/K)^(1/3);
                us0 = (Mad*C/(pi*rho*gamma*Dp^3*(3/4*fD*(1+.15*(gamma/2*(Dp*us00/nu)^2)^.678)+fM )))^(1/2);
                fu = @(us) us - (Mad*C/(pi*rho*gamma*Dp^3*(3/4*fD*(1+.15*(gamma/2*(Dp*us/nu)^2)^.678)+fM )))^(1/2);
                options = optimset('TolFun',1e-15,'Tolx',1e-15);
                ustar = fzero(fu,us0,options);
            case 'Max moment (Babak)' % Mine <<<<<<<<<<<Max moment & nonlinear drag
                gamma = 1.72; %burst
                fD = 1.7009;
                fM = .943993;
                Mmax = 2.707*(Wa^4*Dp^5/K)^(1/3);
                us0 = (Mmax*C/(pi*rho*gamma*Dp^3*(3/4*fD*(1+.15*(gamma/2*(Dp*us00/nu)^2)^.678)+fM )))^(1/2);
                fu = @(us) us - (Mmax*C/(pi*rho*gamma*Dp^3*(3/4*fD*(1+.15*(gamma/2*(Dp*us/nu)^2)^.678)+fM )))^(1/2);
                options = optimset('TolFun',1e-15,'Tolx',1e-15);
                ustar = fzero(fu,us0,options);
        end %switch
%************* Burst/Inrush - Sliding ***********
    case 'TPL-Sliding-Burst-Smooth'
        V = (25*A^2*Dp/(288*z0^7*K^2))^(1/3);
        switch mode
            case 'No gravitational force' %??????
                us0 = ((.5*exp(.124*(V-.01)^.439) + .2*V)*.55*C*mus*Wa/(rho*Dp)/(2.43 + .07*Dp*us00/nu))^.5;
                fu = @(us) us - ((.5*exp(.124*(V-.01)^.439) + .2*V)*.55*C*mus*Wa/(rho*Dp)/(2.43 + .07*Dp*us/nu))^.5;
                options = optimset('TolFun',1e-10,'Tolx',1e-10);
                ustar = fzero(fu,us0,options);
        end %switch
    case 'JKR-Sliding-Burst-Smooth'
        switch mode
            case 'No gravitational force' %??????
                mus
                us0 = .64*(C*mus*Wa/(rho*Dp))^.5 / (2.43 + .07*Dp*us00/nu)^.5
                fu = @(us) us - .64*(C*mus*Wa/(rho*Dp))^.5 / (2.43 + .07*Dp*us/nu)^.5;
                options = optimset('TolFun',1e-10,'Tolx',1e-10);
                [ustar, fval, exitflag, output] = fzero(fu,us0,options);
        end %switch
%*******************************************
%---------------- ROUGH -------------------
%*******************************************
    case 'JKR-Rolling-Sublayer-Rough'
        switch mode
            case 'With lift/gravitational forces'
                beta = parameters.betap * Dp;
                Deltac = parameters.Deltac;
                fpo = theory_pulloff('JKR-asperity',parameters);
                parameters.fpo = fpo;
                deltac = (fpo^2/(3*K^2*beta))^(1/3);
                sigma = deltac/Deltac;
                parameters.sigma = sigma;
                parameters.H0 = 0;
                L = theory_rough_L(0,parameters);
                N = .1/(sigma*beta);
                parameters.N = N;
                a = theory_cradius('JKR-Rough',parameters);
                Fmg = rho_p * 4/3*pi*(Dp/2)^3 * 9.81;
%                 Fmg = 0; %<<<<<
                us0 = (a*C*(pi*a^2*N*fpo*exp(-.6/Deltac^2) + Fmg) / (rho*L*Dp^2*(5.05*pi + 1.95*a*C*us00/nu)))^.5;
                fu = @(us) us - (a*C*(pi*a^2*N*fpo*exp(-.6/Deltac^2) + Fmg) / (rho*L*Dp^2*(5.05*pi + 1.95*a*C*us/nu)))^.5;
                options = optimset('TolFun',1e-10,'Tolx',1e-10);
                ustar = fzero(fu,us0,options);
            case 'Sublayer-(Iman)'
                gamma = 1.1365;
                f = 1.7009;
                fm = .943993;
                Deltac = parameters.Deltac;
                us0 = (.889*Wa^(2/3)*C^.5*exp(-1.2/Deltac^2))/...
                    (K^(1/6)*(gamma*(.75*f*(1+.15*(gamma/2*(us00*Dp/nu)^2)^.678)+fm))^.5*Dp^(2/3));
                fu = @(us) us - (.889*Wa^(2/3)*C^.5*exp(-1.2/Deltac^2))/...
                    (K^(1/6)*(gamma*(.75*f*(1+.15*(gamma/2*(us*Dp/nu)^2)^.678)+fm))^.5*Dp^(2/3));
                options = optimset('TolFun',1e-10,'Tolx',1e-10);
                ustar = fzero(fu,us0,options);
        end %switch
    case 'JKR-Sliding-Sublayer-Rough'
        switch mode
            case 'With lift/gravitational forces'
                beta = parameters.betap * Dp;
                Deltac = parameters.Deltac;
                fpo = theory_pulloff('JKR-asperity',parameters);
                parameters.fpo = fpo;
                deltac = (fpo^2/(3*K^2*beta))^(1/3);
                sigma = deltac/Deltac;
                parameters.sigma = sigma;
                parameters.H0 = 0;
                L = theory_rough_L(0,parameters);
                N = .1/(sigma*beta);
                parameters.N = N;
                a = theory_cradius('JKR-Rough',parameters);
                Fmg = rho_p * 4/3*pi*(Dp/2)^3 * 9.81;
%                 Fmg = 0; %<<<<<
                us0 = (mus*C*(pi*a^2*N*fpo*exp(-.6/Deltac^2) + Fmg) / (rho*L*Dp*(5.8*pi + 1.95*mus*C*Dp*us00/nu)))^.5;
                fu = @(us) us - (mus*C*(pi*a^2*N*fpo*exp(-.6/Deltac^2) + Fmg) / (rho*L*Dp*(5.8*pi + 1.95*mus*C*Dp*us/nu)))^.5;
                options = optimset('TolFun',1e-10,'Tolx',1e-10);
                [ustar, fval, exitflag, output] = fzero(fu,us0,options);
        end %switch
    case 'JKR-Rolling-Burst-Rough'
        switch mode
            case 'With lift/gravitational forces'
                beta = parameters.betap * Dp;
                Deltac = parameters.Deltac;
                fpo = theory_pulloff('JKR-asperity',parameters);
                parameters.fpo = fpo;
                deltac = (fpo^2/(3*K^2*beta))^(1/3);
                sigma = deltac/Deltac;
                parameters.sigma = sigma;
                parameters.H0 = 0;
                L = theory_rough_L(0,parameters);
                N = .1/(sigma*beta);
                parameters.N = N;
                a = theory_cradius('JKR-Rough',parameters);
                Fmg = rho_p * 4/3*pi*(Dp/2)^3 * 9.81;
                Fmg = 0; %<<<<<
                us0 = (a*C*(pi*a^2*N*fpo*exp(-.6/Deltac^2) + Fmg) / (rho*L*Dp^2*(1.72 + .1*us00*L/nu)* (5.04*pi + 1.95*a*C*us00/nu*(1.72 + .2*us00*L/nu)^.5)))^.5;
                fu = @(us) us - (a*C*(pi*a^2*N*fpo*exp(-.6/Deltac^2) + Fmg) / (rho*L*Dp^2*(1.72 + .1*us*L/nu)* (5.04*pi + 1.95*a*C*us/nu*(1.72 + .2*us*L/nu)^.5)))^.5;
                options = optimset('TolFun',1e-10,'Tolx',1e-10);
                ustar = fzero(fu,us0,options);
            case 'Burst-(Iman)'
                if ~isempty(glob_gamma)
                    gamma = glob_gamma;
                else
                    gamma = 1.72;
                end %if
                f = 1.7009;
                fm = .943993;
                Deltac = parameters.Deltac;
                us0 = (.889*Wa^(2/3)*C^.5*exp(-1.2/Deltac^2))/...
                    (K^(1/6)*(gamma*(.75*f*(1+.15*(gamma/2*(us00*Dp/nu)^2)^.678)+fm))^.5*Dp^(2/3));
                fu = @(us) us - (.889*Wa^(2/3)*C^.5*exp(-1.2/Deltac^2))/...
                    (K^(1/6)*(gamma*(.75*f*(1+.15*(gamma/2*(us*Dp/nu)^2)^.678)+fm))^.5*Dp^(2/3));
                options = optimset('TolFun',1e-10,'Tolx',1e-10);
                ustar = fzero(fu,us0,options);
            case 'Burst/Capillary (Babak)'
                if ~isempty(glob_gamma)
                    gamma = glob_gamma;
                else
                    gamma = 1.72;
                end %if
                fD = 1.7009;
                fM = .943993;
                s_tension = .0735; % Water surface tension (N/m)
                Deltac = parameters.Deltac;
                Fc = 2*pi*s_tension*Dp;
                Mad = ((3*pi)^(4/3) / 8) * (Dp^5*Wa^4/K)^(1/3) * exp(-2.4/Deltac^2);
                Mres = Mad + Fc*(Mad*Dp/(2*K))^(1/4);
                
                us0 = (Mres*C/(pi*rho*gamma*Dp^3*(3/4*fD*(1+.15*(gamma/2*(Dp*us00/nu)^2)^.678)+fM )))^(1/2);
                fu = @(us) us - (Mres*C/(pi*rho*gamma*Dp^3*(3/4*fD*(1+.15*(gamma/2*(Dp*us/nu)^2)^.678)+fM )))^(1/2);
                options = optimset('TolFun',1e-15,'Tolx',1e-15);
                ustar = fzero(fu,us0,options);
                
                
        end %switch
    case 'JKR-Sliding-Burst-Rough'
        switch mode
            case 'With lift/gravitational forces'
                beta = parameters.betap * Dp;
                Deltac = parameters.Deltac;
                fpo = theory_pulloff('JKR-asperity',parameters);
                parameters.fpo = fpo;
                deltac = (fpo^2/(3*K^2*beta))^(1/3);
                sigma = deltac/Deltac;
                parameters.sigma = sigma;
                parameters.H0 = 0;
                L = theory_rough_L(0,parameters);
                N = .1/(sigma*beta);
                parameters.N = N;
                a = theory_cradius('JKR-Rough',parameters);
                Fmg = rho_p * 4/3*pi*(Dp/2)^3 * 9.81;
%                 Fmg = 0; %<<<<<
                us0 = (mus*C*(pi*a^2*N*fpo*exp(-.6/Deltac^2) + Fmg) / (rho*Dp*L*(1.72 + .1*us00*L/nu)*(5.8*pi + 1.95*mus*C*Dp*us00/nu*(1.72 + .2*us00*L/nu)^.5)))^.5;
                fu = @(us) us - (mus*C*(pi*a^2*N*fpo*exp(-.6/Deltac^2) + Fmg) / (rho*Dp*L*(1.72 + .1*us*L/nu)*(5.8*pi + 1.95*mus*C*Dp*us/nu*(1.72 + .2*us*L/nu)^.5)))^.5;
                options = optimset('TolFun',1e-10,'Tolx',1e-10);
                [ustar, fval, exitflag, output] = fzero(fu,us0,options);
        end %switch
    case 'JKR-Rolling-Sublayer-Rough-Bumpy'
        switch mode
            case 'No lift/gravitational forces'
                Nbump = parameters.Nbump;
%                 n_b = parameters.n_b;
                n_u = parameters.n_u;
                alpha = parameters.alpha;
                Rb_Dp = parameters.Rb_Dp;
%                 Rb_Dp = 1/(nb*nu*sqrt(Nbump))
                n_b = 1/(n_u*Rb_Dp*sqrt(Nbump));
                parameters.n_b = n_b;
                Rb = Rb_Dp * Dp;
                if (isfield(parameters,'sigma') == 0)
                    Deltac = parameters.Deltac;
                    beta = parameters.betap * Dp;
                    fpo = theory_pulloff('JKR-asperity',parameters);
                    parameters.fpo = fpo;
                    deltac = (fpo^2/(3*K^2*beta))^(1/3);
                    sigma = deltac/Deltac;
                    N = .1/(sigma*beta);
                    glob_sigma(end+1) = sigma;
                    glob_N(end+1) = N;
                    assignin('base','sig',glob_sigma)
                    assignin('base','N',glob_N)
                    parameters.sigma = sigma;
                end %if
%                 N = .1/(sigma*beta);
%                 parameters.N = N;
                Fpo_b = theory_pulloff('JKR-bump',parameters);
                f = 1.7009;
                fm = .943993;
                us0 = (sqrt(3)*Fpo_b*Rb*n_b*C/(1.1365*pi*rho*Dp^3*cos(alpha)*(.75*f*(1+.15*(.56825*(Dp*us00/nu)^2)^.678 )+fm)))^.5;
                fu = @(us) us - (sqrt(3)*Fpo_b*Rb*n_b*C/(1.1365*pi*rho*Dp^3*cos(alpha)* (.75*f*(1+.15*(.56825*(Dp*us/nu)^2)^.678 )+fm)))^.5;
                options = optimset('TolFun',1e-10,'Tolx',1e-10);
                ustar = fzero(fu,us0,options);
            case 'New' % It is not validate yet!!!! (Imans thesis)
                Nbump = parameters.Nbump;
%                 n_b = parameters.n_b;
                n_u = parameters.n_u;
                alpha = parameters.alpha;
%                 Rb_Dp = 1/(n_b*n_u*sqrt(Nbump));
                Rb_Dp = parameters.Rb_Dp;
                n_b = 1/(n_u*Rb_Dp*sqrt(Nbump));
                parameters.n_b = n_b;
                Rb = Rb_Dp * Dp;
%                 if (isfield(parameters,'sigma') == 0)
%                     Deltac = parameters.Deltac;
%                     beta = parameters.betap * Dp;                
%                     fpo = theory_pulloff('JKR-asperity',parameters);
%                     parameters.fpo = fpo;
%                     deltac = (fpo^2/(3*K^2*beta))^(1/3);
%                     sigma = deltac/Deltac;
%                     parameters.sigma = sigma;
%                 end %if
%                 N = .1/(sigma*beta);
%                 parameters.N = N;
                Fpo_b = theory_pulloff('JKR-bump-new',parameters);
                f = 1.7009;
                fm = .943993;
                us0 = (sqrt(3)*Fpo_b*Rb*n_b*C/(1.1365*pi*rho*Dp^3*cos(alpha)*(.75*f*(1+.15*(.56825*(Dp*us00/nu)^2)^.678 )+fm)))^.5;
                fu = @(us) us - (sqrt(3)*Fpo_b*Rb*n_b*C/(1.1365*pi*rho*Dp^3*cos(alpha)* (.75*f*(1+.15*(.56825*(Dp*us/nu)^2)^.678 )+fm)))^.5;
                options = optimset('TolFun',1e-10,'Tolx',1e-10);
                ustar = fzero(fu,us0,options);
        end %switch mode
    case 'JKR-Rolling-Burst-Rough-Bumpy'
        switch mode
            case 'No lift/gravitational forces'
                Nbump = parameters.Nbump;
%                 n_b = parameters.n_b;
                n_u = parameters.n_u;
                alpha = parameters.alpha;
%                 Rb_Dp = 1/(n_b*n_u*sqrt(Nbump));
                Rb_Dp = parameters.Rb_Dp;
                n_b = 1/(n_u*Rb_Dp*sqrt(Nbump));
                parameters.n_b = n_b;
                Rb = Rb_Dp * Dp;
                if (isfield(parameters,'sigma') == 0)
                    Deltac = parameters.Deltac;
                    beta = parameters.betap * Dp;                
                    fpo = theory_pulloff('JKR-asperity',parameters);
                    parameters.fpo = fpo;
                    deltac = (fpo^2/(3*K^2*beta))^(1/3);
                    sigma = deltac/Deltac;
                    parameters.sigma = sigma;
                end %if
%                 N = .1/(sigma*beta);
%                 parameters.N = N;
                Fpo_b = theory_pulloff('JKR-bump',parameters);
                f = 1.7009;
                fm = .943993;
                us0 = (sqrt(3)*Fpo_b*Rb*n_b*C/(1.72*pi*rho*Dp^3*cos(alpha)*(.75*f*(1+.15*(.86*(Dp*us00/nu)^2)^.678)+fm)))^.5;
                fu = @(us) us - (sqrt(3)*Fpo_b*Rb*n_b*C/(1.72*pi*rho*Dp^3*cos(alpha)*(.75*f*(1+.15*(.86*(Dp*us/nu)^2)^.678)+fm)))^.5;
                options = optimset('TolFun',1e-10,'Tolx',1e-10);
                ustar = fzero(fu,us0,options);
            case 'New'
                Nbump = parameters.Nbump;
%                 n_b = parameters.n_b;
                n_u = parameters.n_u;
                alpha = parameters.alpha;
%                 Rb_Dp = 1/(n_b*n_u*sqrt(Nbump));
                Rb_Dp = parameters.Rb_Dp;
                n_b = 1/(n_u*Rb_Dp*sqrt(Nbump));
                parameters.n_b = n_b;
                Rb = Rb_Dp * Dp;
%                 if (isfield(parameters,'sigma') == 0)
%                     Deltac = parameters.Deltac;
%                     beta = parameters.betap * Dp;                
%                     fpo = theory_pulloff('JKR-asperity',parameters);
%                     parameters.fpo = fpo;
%                     deltac = (fpo^2/(3*K^2*beta))^(1/3);
%                     sigma = deltac/Deltac;
%                     parameters.sigma = sigma;
%                 end %if
%                 N = .1/(sigma*beta);
%                 parameters.N = N;
                Fpo_b = theory_pulloff('JKR-bump-new',parameters);
                f = 1.7009;
                fm = .943993;
                us0 = (sqrt(3)*Fpo_b*Rb*n_b*C/(1.72*pi*rho*Dp^3*cos(alpha)*(.75*f*(1+.15*(.86*(Dp*us00/nu)^2)^.678)+fm)))^.5;
                fu = @(us) us - (sqrt(3)*Fpo_b*Rb*n_b*C/(1.72*pi*rho*Dp^3*cos(alpha)*(.75*f*(1+.15*(.86*(Dp*us/nu)^2)^.678)+fm)))^.5;
                options = optimset('TolFun',1e-10,'Tolx',1e-10);
                ustar = fzero(fu,us0,options);
        end %switch mode
        
end %switch model

