function out = SE1P_Rotlet_fourier_space(x, t, opt)
%SE1P_Rotlet_fourier_space  Compute Fourier-space part of the
%1-periodic Ewald sum for the rotlet potential, using the
%Spectral Ewald method.
%
%   out = SE1P_Rotlet_fourier_space(x, t, opt)
%
%   Computes the rotlet velocity field from N point torques.
%   Periodicity is assumed in the first direction (x), while the
%   two remaining directions (y and z) are assumed free.
%
%   Input parameters:
%       :param x: source locations (N×3)
%       :param t: source torques (N×3)
%       :param opt: option struct, see "help SE1P_base_fourier_space"
%
%   :returns: **out** -- output struct, see "help SE1P_base_fourier_space"

opt.kernel = 'rotlet';
if ~isfield(opt,'source_quantity'), opt.source_quantity = sum(t(:).^2); end
out = SE1P_base_fourier_space(x, t, opt);

end
