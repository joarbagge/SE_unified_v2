function out = SE1P_Laplace_fourier_space(x, f, opt)
%SE1P_Laplace_fourier_space  Compute Fourier-space part of the
%1-periodic Ewald sum for the electrostatic potential, using the
%Spectral Ewald method.
%
%   out = SE1P_Laplace_fourier_space(x, f, opt)
%
%   Computes the electrostatic potential from N point charges.
%   Periodicity is assumed in the first direction (x), while the
%   two remaining directions (y and z) are assumed free.
%
%   Input parameters:
%       :param x: source locations (N×3)
%       :param f: source charges (N×1)
%       :param opt: option struct, see "help SE1P_base_fourier_space"
%
%   :returns: **out** -- output struct, see "help SE1P_base_fourier_space"

opt.kernel = 'laplace';
if ~isfield(opt,'source_quantity'), opt.source_quantity = sum(f(:).^2); end
out = SE1P_base_fourier_space(x, f, opt);

end
