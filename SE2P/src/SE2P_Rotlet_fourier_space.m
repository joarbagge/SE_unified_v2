function out = SE2P_Rotlet_fourier_space(x, t, opt)
%SE2P_Rotlet_fourier_space  Compute Fourier-space part of the
%2-periodic Ewald sum for the rotlet potential, using the
%Spectral Ewald method.
%
%   out = SE2P_Rotlet_fourier_space(x, t, opt)
%
%   Computes the rotlet velocity field from N point torques.
%   Periodicity is assumed in the first two directions (x and y),
%   while the third direction (z) is assumed free.
%
%   Input parameters:
%       :param x: source locations (N×3)
%       :param t: source torques (N×3)
%       :param opt: option struct, see "help SE2P_base_fourier_space"
%
%   :returns: **out** -- output struct, see "help SE2P_base_fourier_space"

opt.kernel = 'rotlet';
if ~isfield(opt,'source_quantity'), opt.source_quantity = sum(t(:).^2); end
out = SE2P_base_fourier_space(x, t, opt);

end
