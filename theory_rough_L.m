function L = theory_rough_L(cas,parameters)

Dp = parameters.Dp;
switch cas
    case 0
        L = Dp/2;
    case 1
        sigma = parameters.sigma;
        H0 = parameters.H0;
        L = Dp/2 + 2.76*sigma + H0;
end %switch
