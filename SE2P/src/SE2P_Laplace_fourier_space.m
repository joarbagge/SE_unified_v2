function out = SE2P_Laplace_fourier_space(x, f, opt)
%SE2P_Laplace_fourier_space  Compute Fourier-space part of the
%2-periodic Ewald sum for the electrostatic potential, using the
%Spectral Ewald method.
%
%   out = SE2P_Laplace_fourier_space(x, f, opt)
%
%   Computes the electrostatic potential from N point charges.
%   Periodicity is assumed in the first two directions (x and y),
%   while the third direction (z) is assumed free.
%
%   Input parameters:
%       :param x: source locations (N×3)
%       :param f: source charges (N×1)
%       :param opt: option struct, see "help SE2P_base_fourier_space"
%
%   :returns: **out** -- output struct, see "help SE2P_base_fourier_space"

opt.kernel = 'laplace';
if ~isfield(opt,'source_quantity'), opt.source_quantity = sum(f(:).^2); end
out = SE2P_base_fourier_space(x, f, opt);

end
