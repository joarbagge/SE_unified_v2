function u_self = SE_Laplace_self_term(x, f, opt)
%SE_Laplace_self_term  Compute the self term of the Ewald sum for
%the electrostatic potential.
%
%   u_self = SE_Laplace_self_term(x, f, opt)
%
%   Input parameters:
%       :param x: source locations (N×3)
%       :param f: source charges (N×1)
%       :param opt: option struct, see below
%
%   Valid options for the opt struct:
%       :param opt.xi: Ewald decomposition parameter (required)
%
%   :returns: **u_self** -- potential self term

assert(isfield(opt, 'xi'), 'Ewald decomposition parameter "xi" must be given in opt struct');

u_self = (-2*opt.xi/sqrt(pi)) * f;

end
