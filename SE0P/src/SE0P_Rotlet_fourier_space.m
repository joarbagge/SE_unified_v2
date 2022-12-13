function out = SE0P_Rotlet_fourier_space(x, t, opt, pre_kernel)
%SE0P_Rotlet_fourier_space  Compute Fourier-space part of the
%0-periodic (free-space) Ewald sum for the rotlet potential,
%using the Spectral Ewald method.
%
%   out = SE0P_Rotlet_fourier_space(x, t, opt)
%   out = SE0P_Rotlet_fourier_space(x, t, opt, pre_kernel)
%
%   Computes the rotlet velocity field from N point torques.
%
%   Input parameters:
%       :param x: source locations (N×3)
%       :param t: source torques (N×3)
%       :param opt: option struct, see "help SE0P_base_fourier_space"
%       :param pre_kernel: optional, output from SE0P_Rotlet_precompute_kernel_fft
%
%   :returns: **out** -- output struct, see "help SE0P_base_fourier_space"

if nargin < 4 || isempty(pre_kernel)
  pre_kernel = SE0P_Rotlet_precompute_kernel_fft(opt);
end

opt.kernel = 'rotlet';
out = SE0P_base_fourier_space(x, t, opt, pre_kernel);

end
