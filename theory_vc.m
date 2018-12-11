function ustar = theory_vc(model,mode,parameters)
ustar = -1;
k = 1;
pass = 0;
% ustep = .0001;
% us = -ustep;
tolerance = 1e-15;
us = 1;
ustep = 1;
% if ~isfield(parameters,'assumption')
%     parameters.assumption.nldrag = 0;
%     parameters.assumption.liftforce = 0;
%     parameters.assumption.gravityforce = 0;
% else
%     if ~isfield(parameters.assumption,'nldrag')
%         parameters.assumption.nldrag = 0;
%     end %if
%     if ~isfield(parameters.assumption,'liftforce')
%         parameters.assumption.liftforce = 0;
%     end %if
%     if ~isfield(parameters.assumption,'gravityforce')
%         parameters.assumption.gravityforce = 0;
%     end %if
% end %if
switch model
% -------------------------------------------------------------------------
% ******************************** Rolling ********************************
% -------------------------------------------------------------------------
% *********** Rolling - Smooth - Sublayer ***************
    case 'TPL-Rolling-Sublayer-Smooth'
        switch mode
            case 'Detached moment'
                parameters.a = theory_cradius('TPL',parameters);
                parameters.Fpo = theory_pulloff('TPL',parameters);  
                Mr = theory_rmoment('Detached moment-Smooth',parameters);
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    parameters.Fd = theory_drag('Stokes-Sublayer-Smooth',parameters);
                    parameters.Fl = theory_lift('Saffman-Sublayer-Smooth',parameters);
                    parameters.M = theory_momentp('Sublayer-Smooth',parameters);
                    Mh = theory_hmoment('Smooth',parameters);
                    if ((Mh > Mr) && (ustep > 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if ((Mh < Mr) && (ustep < 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if (abs((Mr-Mh)/Mr) <= tolerance)
                        ustar = us;
                        k = 0;
                        return
                    end %if                    
%                     if (Mh >= Mr)
%                         k = 0;
%                         ustar = us;
%                         return
%                     end %if
                end %while
        end %switch mode
    case 'JKR-Rolling-Sublayer-Smooth'
        switch mode
            case 'Detached moment'
                parameters.a = theory_cradius('JKR',parameters);
                parameters.Fpo = theory_pulloff('JKR',parameters);
                Mr = theory_rmoment('Detached moment-Smooth',parameters);
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    parameters.Fd = theory_drag('Stokes-Sublayer-Smooth',parameters);
                    if (parameters.assumption.nldrag == 1)
                        parameters.Fd = theory_drag('Nonlinear-Sublayer',parameters);
                    end %if
                    if (parameters.assumption.liftforce == 1)
                        parameters.Fl = theory_lift('Saffman-Sublayer-Smooth',parameters);
                    else
                        parameters.Fl = 0;
                    end %if
                    parameters.M = theory_momentp('Sublayer-Smooth',parameters);
                    Mh = theory_hmoment('Smooth',parameters);
                    if ((Mh > Mr) && (ustep > 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if ((Mh < Mr) && (ustep < 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if (abs((Mr-Mh)/Mr) <= tolerance)
                        ustar = us;
                        k = 0;
                        return
                    end %if
%                     if (Mh >= Mr)
%                         k = 0;
%                         ustar = us;
%                         return
%                     end %if
                end %while
            case 'Max moment'
                parameters.a = theory_cradius('JKR-Mmax',parameters);
                Mr = theory_rmoment('MmaxJKR-Smooth',parameters);
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    parameters.Fd = theory_drag('Stokes-Sublayer-Smooth',parameters);
                    if (parameters.assumption.nldrag == 1)
                        parameters.Fd = theory_drag('Nonlinear-Sublayer',parameters);
                    end %if
                    if (parameters.assumption.liftforce == 1)
                        parameters.Fl = theory_lift('Saffman-Sublayer-Smooth',parameters);
                    else
                        parameters.Fl = 0;
                    end %if
                    parameters.M = theory_momentp('Sublayer-Smooth',parameters);
                    Mh = theory_hmoment('Smooth',parameters);
                    if ((Mh > Mr) && (ustep > 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if ((Mh < Mr) && (ustep < 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if (abs((Mr-Mh)/Mr) <= tolerance)
                        ustar = us;
                        k = 0;
                        return
                    end %if                         
%                     if (Mh >= Mr)
%                         k = 0;
%                         ustar = us;
%                         return
%                     end %if
                end %while
        end %switch mode
    case 'DMT-Rolling-Sublayer-Smooth'
        switch mode
            case 'Max moment'
                parameters.a = theory_cradius('DMT-Mmax',parameters);
                Mr = theory_rmoment('MmaxDMT-Smooth',parameters);
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    parameters.Fd = theory_drag('Stokes-Sublayer-Smooth',parameters);
                    if (parameters.assumption.nldrag == 1)
                        parameters.Fd = theory_drag('Nonlinear-Sublayer',parameters);
                    end %if
                    if (parameters.assumption.liftforce == 1)
                        parameters.Fl = theory_lift('Saffman-Sublayer-Smooth',parameters);
                    else
                        parameters.Fl = 0;
                    end %if
                    parameters.M = theory_momentp('Sublayer-Smooth',parameters);
                    Mh = theory_hmoment('Smooth',parameters);
                    if ((Mh > Mr) && (ustep > 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if ((Mh < Mr) && (ustep < 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if (abs((Mr-Mh)/Mr) <= tolerance)
                        ustar = us;
                        k = 0;
                        return
                    end %if
%                     if (Mh >= Mr)
%                         k = 0;
%                         ustar = us;
%                         return
%                     end %if
                end %while
        end %switch
%*********** Rolling - Smooth - Burst/Inrush ***********
    case 'TPL-Rolling-Burst-Smooth'
        parameters.a = theory_cradius('TPL',parameters);
        parameters.Fpo = theory_pulloff('TPL',parameters);
        switch mode
            case 'Detached moment'
                Mr = theory_rmoment('Detached moment-Smooth',parameters);               
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    parameters.Fd = theory_drag('Stokes-Burst-Smooth',parameters);
                    if (parameters.assumption.nldrag == 1)
                        parameters.Fd = theory_drag('Nonlinear-Burst',parameters);
                    end %if
                    parameters.Fl = theory_lift('Saffman-Burst-Smooth',parameters);
                    parameters.M = theory_momentp('Burst-Smooth',parameters);
                    Mh = theory_hmoment('Smooth',parameters);
                    if ((Mh > Mr) && (ustep > 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if ((Mh < Mr) && (ustep < 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if (abs((Mr-Mh)/Mr) <= tolerance)
                        ustar = us;
                        k = 0;
                        return
                    end %if
%                     if (Mh >= Mr)
%                         k = 0;
%                         ustar = us;
%                         return
%                     end %if
                end %while
        end %switch
    case 'JKR-Rolling-Burst-Smooth'
        switch mode
            case 'Detached moment'
                parameters.a = theory_cradius('JKR',parameters);
                parameters.Fpo = theory_pulloff('JKR',parameters);
                Mr = theory_rmoment('Detached moment-Smooth',parameters);
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    parameters.Fd = theory_drag('Stokes-Burst-Smooth',parameters);
                    if (parameters.assumption.nldrag == 1)
%                         disp ('nonlinear')
                        parameters.Fd = theory_drag('Nonlinear-Burst',parameters);
                    end %if
                    
                    if (parameters.assumption.liftforce == 1)
                        parameters.Fl = theory_lift('Saffman-Burst-Smooth',parameters);
                    else
                        parameters.Fl = 0;
                    end %if
                    parameters.M = theory_momentp('Burst-Smooth',parameters);
                    Mh = theory_hmoment('Smooth',parameters);
                    if ((Mh > Mr) && (ustep > 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if ((Mh < Mr) && (ustep < 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if (abs((Mr-Mh)/Mr) <= tolerance)
                        ustar = us;
                        k = 0;
                        return
                    end %if
%                     if (Mh >= Mr)
%                         k = 0;
%                         ustar = us;
%                         return
%                     end %if
                end %while
            case 'Max moment'
                parameters.a = theory_cradius('JKR-Mmax',parameters);
                Mr = theory_rmoment('MmaxJKR-Smooth',parameters);
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    parameters.Fd = theory_drag('Stokes-Burst-Smooth',parameters);
                    if (parameters.assumption.nldrag == 1)
                        parameters.Fd = theory_drag('Nonlinear-Burst',parameters);
                    end %if
                    if (parameters.assumption.liftforce == 1)
                        parameters.Fl = theory_lift('Saffman-Burst-Smooth',parameters);
                    else
                        parameters.Fl = 0;
                    end %if
                    parameters.M = theory_momentp('Burst-Smooth',parameters);
                    Mh = theory_hmoment('Smooth',parameters);
                    if ((Mh > Mr) && (ustep > 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if ((Mh < Mr) && (ustep < 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if (abs((Mr-Mh)/Mr) <= tolerance)
                        ustar = us;
                        k = 0;
                        return
                    end %if
%                     if (Mh >= Mr)
%                         k = 0;
%                         ustar = us;
%                         return
%                     end %if
                end %while
        end %switch
    case 'DMT-Rolling-Burst-Smooth'
        switch mode
            case 'Max moment'
                parameters.a = theory_cradius('DMT-Mmax',parameters);
                Mr = theory_rmoment('MmaxDMT-Smooth',parameters);
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    parameters.Fd = theory_drag('Stokes-Burst-Smooth',parameters);
                    if (parameters.assumption.nldrag == 1)
                        parameters.Fd = theory_drag('Nonlinear-Burst',parameters);
                    end %if
                    if (parameters.assumption.liftforce == 1)
                        parameters.Fl = theory_lift('Saffman-Burst-Smooth',parameters);
                    else
                        parameters.Fl = 0;
                    end %if
                    parameters.M = theory_momentp('Burst-Smooth',parameters);
                    Mh = theory_hmoment('Smooth',parameters);
                    if ((Mh > Mr) && (ustep > 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if ((Mh < Mr) && (ustep < 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if (abs((Mr-Mh)/Mr) <= tolerance)
                        ustar = us;
                        k = 0;
                        return
                    end %if
%                     if (Mh >= Mr)
%                         k = 0;
%                         ustar = us;
%                         return
%                     end %if
                end %while
        end %switch
% ______________________________________
% ------------ Rolling - Rough --------
% ______________________________________
%********* Rolling - Rough - Sublayer ********
    case 'JKR-Rolling-Sublayer-Rough'
        Dp = parameters.Dp;
        beta = parameters.betap * Dp;
        Deltac = parameters.Deltac;
        switch mode
            case 'Detached moment'
                parameters.fpo = theory_pulloff('JKR-asperity',parameters);
                deltac = ((parameters.fpo)^2/(3*(parameters.K)^2*beta))^(1/3);
                sigma = deltac/Deltac;
                parameters.sigma = sigma;
                parameters.H0 = 0;
                parameters.L = theory_rough_L(0,parameters);
                N = .1/(sigma*beta);
                parameters.N = N;                
                parameters.a = theory_cradius('JKR-Rough',parameters);
                parameters.FM = theory_pulloff('JKR-Rough',parameters);
                Mr = theory_rmoment('Detached moment-Rough',parameters);
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    parameters.Fd = theory_drag('Stokes-Sublayer-Rough',parameters);
                    if (parameters.assumption.nldrag == 1)
                        parameters.Fd = theory_drag('Nonlinear-Sublayer',parameters);
                    end %if
                    parameters.Fl = theory_lift('Saffman-Sublayer-Rough',parameters);
                    parameters.M = theory_momentp('Sublayer-Rough',parameters);
                    Mh = theory_hmoment('Rough',parameters);
                    if ((Mh > Mr) && (ustep > 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if ((Mh < Mr) && (ustep < 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if (abs((Mr-Mh)/Mr) <= tolerance)
                        ustar = us;
                        k = 0;
                        return
                    end %if
%                     if (Mh >= Mr)
%                         k = 0;
%                         ustar = us;
%                         return
%                     end %if
                end %while
        end %switch
%*********** Rolling - Rough - Burst ************
    case 'JKR-Rolling-Burst-Rough'
        Dp = parameters.Dp;
        beta = parameters.betap * Dp;
        Deltac = parameters.Deltac;
        switch mode
            case 'Detached moment'
                parameters.fpo = theory_pulloff('JKR-asperity',parameters);
                deltac = ((parameters.fpo)^2/(3*(parameters.K)^2*beta))^(1/3);
                sigma = deltac/Deltac;
                parameters.sigma = sigma;
                parameters.H0 = 0;
                parameters.L = theory_rough_L(0,parameters);
                N = .1/(sigma*beta);
                parameters.N = N;
                parameters.a = theory_cradius('JKR-Rough',parameters);
                parameters.FM = theory_pulloff('JKR-Rough',parameters);
                Mr = theory_rmoment('Detached moment-Rough',parameters);
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    parameters.Fd = theory_drag('Stokes-Burst-Rough',parameters);
                    if (parameters.assumption.nldrag == 1)
                        parameters.Fd = theory_drag('Nonlinear-Burst',parameters);
                    end %if
                    parameters.Fl = theory_lift('Saffman-Burst-Rough',parameters);
                    parameters.M = theory_momentp('Burst-Rough',parameters);
                    Mh = theory_hmoment('Rough',parameters);
                    if ((Mh > Mr) && (ustep > 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if ((Mh < Mr) && (ustep < 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if (abs((Mr-Mh)/Mr) <= tolerance)
                        ustar = us;
                        k = 0;
                        return
                    end %if
%                     if (Mh >= Mr)
%                         k = 0;
%                         ustar = us;
%                         return
%                     end %if
                end %while
            case 'Burst-(Iman)'
                parameters.L = theory_rough_L(0,parameters);
                parameters.a = theory_cradius('JKR-Rough-Iman',parameters);
                parameters.FM = theory_pulloff('JKR-Rough-Iman',parameters);
                Mr = theory_rmoment('Detached moment-Rough',parameters);
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    parameters.Fd = theory_drag('Stokes-Burst-Rough',parameters);
                    if (parameters.assumption.nldrag == 1)
                        parameters.Fd = theory_drag('Nonlinear-Burst',parameters);
                    end %if
                    parameters.Fl = theory_lift('Saffman-Burst-Rough',parameters);
                    parameters.M = theory_momentp('Burst-Rough',parameters);
                    Mh = theory_hmoment('Rough',parameters);
                    if ((Mh > Mr) && (ustep > 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if ((Mh < Mr) && (ustep < 0))
                        ustep = -ustep / 2;
                        pass = pass + 1;
                    end %if
                    if (abs((Mr-Mh)/Mr) <= tolerance)
                        ustar = us;
                        k = 0;
                        return
                    end %if
                end %while
        end %switch
% -------------------------------------------------------------------------
% ******************************** Sliding ********************************
% -------------------------------------------------------------------------
% ______________________________________
% -------------Sliding - Smooth --------
% ______________________________________
%********* Sliding - Smooth - Sublayer **********        
    case 'TPL-Sliding-Sublayer-Smooth'
        parameters.Fpo = theory_pulloff('TPL',parameters);
        mus = parameters.mus;
        switch mode
            case 'With lift/gravitational forces'
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    Fd = theory_drag('Stokes-Sublayer-Smooth',parameters);
                    Fl = theory_lift('Saffman-Sublayer-Smooth',parameters);
                    parameters.Fl = Fl;
                    Fn = theory_normalforce('With lift/gravitational forces',parameters);
                    if (Fd >= mus*Fn)
                        k = 0;
                        ustar = us;
                        return
                    end %if
                end %while
        end %switch
    case 'JKR-Sliding-Sublayer-Smooth'
        parameters.Fpo = theory_pulloff('JKR',parameters);
        mus = parameters.mus;
        switch mode
            case 'With lift/gravitational forces'
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    Fd = theory_drag('Stokes-Sublayer-Smooth',parameters);
                    Fl = theory_lift('Saffman-Sublayer-Smooth',parameters);
                    parameters.Fl = Fl;
                    Fn = theory_normalforce('With lift/gravitational forces',parameters);
                    if (Fd >= mus*Fn)
                        k = 0;
                        ustar = us;
                        return
                    end %if
                end %while
            case 'No gravitational force'
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    Fd = theory_drag('Stokes-Sublayer-Smooth',parameters);
                    Fl = theory_lift('Saffman-Sublayer-Smooth',parameters);
                    parameters.Fl = Fl;
                    Fn = theory_normalforce('No gravitational force',parameters);
                    if (Fd >= mus*Fn)
                        k = 0;
                        ustar = us;
                        return
                    end %if
                end %while
            case 'No lift/gravitational forces'
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    Fd = theory_drag('Stokes-Sublayer-Smooth',parameters);
                    Fl = theory_lift('Saffman-Sublayer-Smooth',parameters);
                    parameters.Fl = Fl;
                    Fn = theory_normalforce('No lift/gravitational forces',parameters);
                    if (Fd >= mus*Fn)
                        k = 0;
                        ustar = us;
                        return
                    end %if
                end %while
        end %switch
%************ Sliding - Smooth - Burst/Inrush **********
    case 'TPL-Sliding-Burst-Smooth'
        parameters.Fpo = theory_pulloff('TPL',parameters);
        mus = parameters.mus;
        switch mode
            case 'With lift/gravitational forces'
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    Fd = theory_drag('Stokes-Burst-Smooth',parameters);
                    Fl = theory_lift('Saffman-Burst-Smooth',parameters);
                    parameters.Fl = Fl;
                    Fn = theory_normalforce('With lift/gravitational forces',parameters);
                    if (Fd >= mus*Fn)
                        k = 0;
                        ustar = us;
                        return
                    end %if
                end %while
        end %switch
    case 'JKR-Sliding-Burst-Smooth'
        parameters.Fpo = theory_pulloff('JKR',parameters);
        mus = parameters.mus;
        switch mode
            case 'With lift/gravitational forces'
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    Fd = theory_drag('Stokes-Burst-Smooth',parameters);
                    Fl = theory_lift('Saffman-Burst-Smooth',parameters);
                    parameters.Fl = Fl;
                    Fn = theory_normalforce('With lift/gravitational forces',parameters);
                    if (Fd >= mus*Fn)
                        k = 0;
                        ustar = us;
                        return
                    end %if
                end %while
            case 'No gravitational force'
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    Fd = theory_drag('Stokes-Burst-Smooth',parameters);
                    Fl = theory_lift('Saffman-Burst-Smooth',parameters);
                    parameters.Fl = Fl;
                    Fn = theory_normalforce('No gravitational force',parameters);
                    if (Fd >= mus*Fn)
                        k = 0;
                        ustar = us;
                        return
                    end %if
                end %while
            case 'No lift/gravitational forces'
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    Fd = theory_drag('Stokes-Burst-Smooth',parameters);
                    Fl = theory_lift('Saffman-Burst-Smooth',parameters);
                    parameters.Fl = Fl;
                    Fn = theory_normalforce('No lift/gravitational forces',parameters);
                    if (Fd >= mus*Fn)
                        k = 0;
                        ustar = us;
                        return
                    end %if
                end %while
        end %switch
% ______________________________________
% ----------- Sliding - Rough ----------
% ______________________________________
%********* Sliding - Rough - Sublayer *********
    case 'JKR-Sliding-Sublayer-Rough'
        Dp = parameters.Dp;
        beta = parameters.betap * Dp;
        Deltac = parameters.Deltac;
        parameters.fpo = theory_pulloff('JKR-asperity',parameters);
        deltac = ((parameters.fpo)^2/(3*(parameters.K)^2*beta))^(1/3);
        sigma = deltac/Deltac;
        parameters.sigma = sigma;
        parameters.H0 = 0;
        parameters.L = theory_rough_L(0,parameters);
        N = .1/(sigma*beta);
        parameters.N = N;                
        parameters.a = theory_cradius('JKR-Rough',parameters);
        parameters.FM = theory_pulloff('JKR-Rough',parameters);
        parameters.Fpo = parameters.FM;
        mus = parameters.mus;
        switch mode
            case 'With lift/gravitational forces'
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    Fd = theory_drag('Stokes-Sublayer-Rough',parameters);
                    Fl = theory_lift('Saffman-Sublayer-Rough',parameters);
                    parameters.Fl = Fl;
                    Fn = theory_normalforce('With lift/gravitational forces',parameters);
                    if (Fd >= mus*Fn)
                        k = 0;
                        ustar = us;
                        return
                    end %if
                end %while
        end %switch
%********* Sliding - Rough - Burst ***************
        case 'JKR-Sliding-Burst-Rough'
        Dp = parameters.Dp;
        beta = parameters.betap * Dp;
        Deltac = parameters.Deltac;
        parameters.fpo = theory_pulloff('JKR-asperity',parameters);
        deltac = ((parameters.fpo)^2/(3*(parameters.K)^2*beta))^(1/3);
        sigma = deltac/Deltac;
        parameters.sigma = sigma;
        parameters.H0 = 0;
        parameters.L = theory_rough_L(0,parameters);
        N = .1/(sigma*beta);
        parameters.N = N;                
        parameters.a = theory_cradius('JKR-Rough',parameters);
        parameters.FM = theory_pulloff('JKR-Rough',parameters);
        parameters.Fpo = parameters.FM;
        mus = parameters.mus;
        switch mode
            case 'With lift/gravitational forces'
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    Fd = theory_drag('Stokes-Burst-Rough',parameters);
                    Fl = theory_lift('Saffman-Burst-Rough',parameters);
                    parameters.Fl = Fl;
                    Fn = theory_normalforce('With lift/gravitational forces',parameters);
                    if (Fd >= mus*Fn)
                        k = 0;
                        ustar = us;
                        return
                    end %if
                end %while
        end %switch
% -------------------------------------------------------------------------
% ******************************** Lifting ********************************
% -------------------------------------------------------------------------
%********** Lifting - Smooth - Sublayer ********
    case 'TPL-Lifting-Sublayer-Smooth'
        parameters.Fpo = theory_pulloff('TPL',parameters);
        switch mode
            case 'With gravitational force'
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    Fl = theory_lift('Saffman-Sublayer-Smooth',parameters);
                    parameters.Fl = Fl;
                    Fn = theory_normalforce('With lift/gravitational forces',parameters);
                    if (Fn <= 0)
                        k = 0;
                        ustar = us;
                        return
                    end %if
                end %while
        end %switch
    case 'JKR-Lifting-Sublayer-Smooth'
        parameters.Fpo = theory_pulloff('JKR',parameters);
        switch mode
            case 'With gravitational force'
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    Fl = theory_lift('Saffman-Sublayer-Smooth',parameters);
                    parameters.Fl = Fl;
                    Fn = theory_normalforce('With lift/gravitational forces',parameters);
                    if (Fn <= 0)
                        k = 0;
                        ustar = us;
                        return
                    end %if
                end %while
        end %switch
%*********** Lifting - Smooth - Burst/Inrush *********
    case 'TPL-Lifting-Burst-Smooth'
        parameters.Fpo = theory_pulloff('TPL',parameters);
        switch mode
            case 'With gravitational force'
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    Fl = theory_lift('Saffman-Burst-Smooth',parameters);
                    parameters.Fl = Fl;
                    Fn = theory_normalforce('With lift/gravitational forces',parameters);
                    if (Fn <= 0)
                        k = 0;
                        ustar = us;
                        return
                    end %if
                end %while
        end %switch
    case 'JKR-Lifting-Burst-Smooth'
        parameters.Fpo = theory_pulloff('JKR',parameters);
        switch mode
            case 'With gravitational force'
                while (k == 1)
                    us = us + ustep;
                    parameters.ustar = us;
                    Fl = theory_lift('Saffman-Burst-Smooth',parameters);
                    parameters.Fl = Fl;
                    Fn = theory_normalforce('With lift/gravitational forces',parameters);
                    if (Fn <= 0)
                        k = 0;
                        ustar = us;
                        return
                    end %if
                end %while
        end %switch
%********** Lifting - Rough - Burst ************

end %switch model

