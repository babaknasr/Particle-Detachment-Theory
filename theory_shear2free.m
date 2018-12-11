% convert shear velocity to free-stream (bulk) veloicty
% uf = uf = theory_shear2free(ref,parameters)
% ref: 
%     Blasius(1913)
%     ...
% parameters. 
%           ustar: shear velocity
%           nu:
%           Dh

function uf = theory_shear2free(ref,parameters)

us = parameters.ustar;
switch ref
    case {'Do not Convert','Do not convert','do not convert','None','none'}
        uf = us;
    case 'Schlichting(1979)'
        uf = (us - .0387)/.0375;
    case 'Davies(1972)'
        nu = parameters.nu;
        Dh = parameters.Dh;
        uf = ((5*us).^8*Dh/nu)^(1/7);
    case 'Blasius(1913)'
        nu = parameters.nu;
        Dh = parameters.Dh;
        uf = (sqrt(2/.0791).*us).^(8/7) * (Dh/nu)^(1/7);
end %switch
