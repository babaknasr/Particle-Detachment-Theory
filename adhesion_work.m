% Work of adhesion from Hamaker Constant
function Wa = adhesion_work(varargin)

A = varargin{1};
if (length(varargin) == 2)
    z0 = varargin{2};
else
    z0 = 4e-10;
end %if

Wa = A/(12*pi*z0^2);