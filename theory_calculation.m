function [Dp, theory_result,cmodels,cmodes] = theory_calculation(theory)
global theory_isprogressbar
if isempty(theory_isprogressbar)
    isprogressbar = 1;
else
    isprogressbar = theory_isprogressbar;
end %if
% Air Properties
par.rho = theory.const.rho;
if isfield(theory.const,'nu')
    par.nu = theory.const.nu;
else
    par.mu = theory.const.mu;
    par.nu = par.mu/par.rho;
end
lambda = theory.const.lambda;
% check the existence of z0
if ~isfield(theory.const,'z0')
    theory.const.z0 = 4e-10;
end %if
% Check the particle density
if isnan(theory.const.rho_p)
    theory.const.rho_p = 0;
end %if
% E1 = theory.const.E1;
% E2 = theory.const.E2;
% nu1 = theory.const.nu1;
% nu2 = theory.const.nu2;
% par.K = comp_modulus(E1,nu1,E2,nu2);
par.K = theory.const.K;
par.Wa = theory.const.Wa;
par.rho_p = theory.const.rho_p;
par.A = theory.const.A;
par.z0 = theory.const.z0;
par.mus = theory.const.mus;
%*****
allmodes = {'Sublayer-(Iman)','Burst-(Iman)','Detached moment','Max moment','With lift/gravitational forces',...
    'With gravitational force','No lift/gravitational forces','No gravitational force','New','Max moment (Babak)','Detached moment (Babak)','Burst/Capillary (Babak)'};
jj = length(allmodes);
j = 1;
models = {};
modes = {};
if (theory.smooth == 1)
    if (theory.sublayer == 1)
        if (theory.jkr == 1)
            if (theory.roll == 1)
                models(j:jj+j-1) = {'JKR-Rolling-Sublayer-Smooth'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %if
            if (theory.slide == 1)
                models(j:jj+j-1) = {'JKR-Sliding-Sublayer-Smooth'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %if
            if (theory.lift == 1)
                models(j:jj+j-1) = {'JKR-Lifting-Sublayer-Smooth'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %if
        end %if
        if (theory.tpl == 1)
            if (theory.roll == 1)
                models(j:jj+j-1) = {'TPL-Rolling-Sublayer-Smooth'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %if
            if (theory.slide == 1)
                models(j:jj+j-1) = {'TPL-Sliding-Sublayer-Smooth'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %if
            if (theory.lift == 1)
                models(j:jj+j-1) = {'TPL-Lifting-Sublayer-Smooth'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %if
        end %if
        if (theory.dmt == 1)
            if (theory.roll == 1)
                models(j:jj+j-1) = {'DMT-Rolling-Sublayer-Smooth'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %if
        end %if
    end %if sublayer
    if (theory.burst == 1)
        if (theory.jkr == 1)
            if (theory.roll == 1)
                models(j:jj+j-1) = {'JKR-Rolling-Burst-Smooth'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %if
            if (theory.slide == 1)
                models(j:jj+j-1) = {'JKR-Sliding-Burst-Smooth'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %if
            if (theory.lift == 1)
                models(j:jj+j-1) = {'JKR-Lifting-Burst-Smooth'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %if
        end %if
        if (theory.tpl == 1)
            if (theory.roll == 1)
                models(j:jj+j-1) = {'TPL-Rolling-Burst-Smooth'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %if
            if (theory.slide == 1)
                models(j:jj+j-1) = {'TPL-Sliding-Burst-Smooth'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %if
            if (theory.lift == 1)
                models(j:jj+j-1) = {'TPL-Lifting-Burst-Smooth'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %if
        end %if
        if (theory.dmt == 1)
            if (theory.roll == 1)
                models(j:jj+j-1) = {'DMT-Rolling-Burst-Smooth'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %if
        end %if
    end %if burst
end %if smooth
% >>>>>>>>>>ROUGH<<<<<<<<<<
if (theory.rough == 1)
    par.betap = theory.roughness.betap;
    par.Deltac = theory.roughness.Deltac;
    if (theory.sublayer == 1)
        if (theory.jkr == 1)
            if (theory.roll == 1)
                models(j:jj+j-1) = {'JKR-Rolling-Sublayer-Rough'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %if roll
            if (theory.slide == 1)
                models(j:jj+j-1) = {'JKR-Sliding-Sublayer-Rough'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %if slide
            if (theory.lift == 1)
                models(j:jj+j-1) = {'JKR-Lifting-Sublayer-Rough'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %if lift
        end %if jkr
    end %if sublayer
    if (theory.burst == 1)
        if (theory.jkr == 1)
            if (theory.roll == 1)
                models(j:jj+j-1) = {'JKR-Rolling-Burst-Rough'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %if roll
            if (theory.slide == 1)
                models(j:jj+j-1) = {'JKR-Sliding-Burst-Rough'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %if slide
            if (theory.lift == 1)
                models(j:jj+j-1) = {'JKR-Lifting-Burst-Rough'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %if lift
        end %if jkr
    end %if burst
end %if rough
% <<<<<<<<< ROUGH/BUMPY >>>>>>>>>>>>
if (theory.rough_bumpy == 1)
    if (isfield(theory.roughness,'sigma') == 0)
        par.betap = theory.roughness.betap;
    else
        par.sigma = theory.roughness.sigma;
    end %if
    par.Deltac = theory.roughness.Deltac;
    par.Nbump = theory.roughness.Nbump;
%     par.n_b = theory.roughness.n_b;
    par.Rb_Dp = theory.roughness.Rb_Dp;
    par.n_u = theory.roughness.n_u;
    par.alpha = theory.roughness.alpha;
    if (theory.sublayer == 1)
        if (theory.jkr == 1)
            if (theory.roll == 1)
                models(j:jj+j-1) = {'JKR-Rolling-Sublayer-Rough-Bumpy'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %roll
        end %jkr
    end %sublayer
    if (theory.burst == 1)
        if (theory.jkr == 1)
            if (theory.roll == 1)
                models(j:jj+j-1) = {'JKR-Rolling-Burst-Rough-Bumpy'};
                modes(j:jj+j-1) = allmodes;
                j = j + jj;
            end %roll
        end %jkr
    end %burst
end %if bumpy

%%% For Function "2"
if ~isfield(theory,'assumption')
    theory.assumption.nldrag = 1;
    theory.assumption.liftforce = 0;
    theory.assumption.gravityforce = 0;
else
    if ~isfield(theory.assumption,'nldrag')
        theory.assumption.nldrag = 1;
    end %if
    if ~isfield(theory.assumption,'liftforce')
        theory.assumption.lift = 0;
    end %if
    if ~isfield(theory.assumption,'gravityforce')
        theory.assumption.gravityforce = 0;
    end %if
end %if
par.assumption = theory.assumption;
%%%% <------

Nsteps = 500;
if (theory.minsize ~= theory.maxsize)
    Dp = linspace(theory.minsize,theory.maxsize,Nsteps);
else
    Dp = theory.minsize;
end %if

theory_result = [];
m = 0;
cmodels = {};
cmodes = {};
if (isprogressbar == 1)
    wb = waitbar(0,'Computing... Please wait...');
end %if
for k=1:length(models)
    ufree = zeros(1,length(Dp));
    for i=1:length(Dp)
        par.C = theory_Cc(Dp(i),lambda,[]);
        par.Dp = Dp(i);
        switch theory.function
            case 1
                ustar = theory_vc_soltani(models{k},modes{k},par);
            case 2
                ustar = theory_vc(models{k},modes{k},par);
        end %switch
        if (ustar == -1)
            break
        end %if
        par2.ustar = ustar;
        par2.nu = par.nu;
        if ((strcmpi(theory.freevelocity,'do not convert') == 0) && (strcmpi(theory.freevelocity,'none') == 0))
            par2.Dh = theory.Dh;
        end %if
        ufree(i) = theory_shear2free(theory.freevelocity,par2);
    end %i
    if (ustar ~= -1)
        m = m + 1;
        cmodels{m} = models{k};
        cmodes{m} = modes{k};
        theory_result(m,:) = ufree;
    end %if
    if (isprogressbar == 1)
        waitbar(k/length(models));
    end %if
end %k
if (isprogressbar == 1)
    delete(wb);
end %if
