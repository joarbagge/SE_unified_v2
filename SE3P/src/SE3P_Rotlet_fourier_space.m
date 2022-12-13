function out = SE3P_Rotlet_fourier_space(x, t, opt)
%SE3P_Rotlet_fourier_space  Compute Fourier-space part of the
%3-periodic Ewald sum for the rotlet potential, using the
%Spectral Ewald method.
%
%   out = SE3P_Rotlet_fourier_space(x, t, opt)
%
%   Computes the rotlet velocity field from N point torques.
%
%   Input parameters:
%       :param x: source locations (N×3)
%       :param t: source torques (N×3)
%       :param opt: option struct, see "help SE3P_base_fourier_space"
%
%   :returns: **out** -- output struct, see "help SE3P_base_fourier_space"

opt.kernel = 'rotlet';
out = SE3P_base_fourier_space(x, t, opt);

end
