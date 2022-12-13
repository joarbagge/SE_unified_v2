function out = SE0P_Laplace_fourier_space(x, f, opt, pre_kernel)
%SE0P_Laplace_fourier_space  Compute Fourier-space part of the
%0-periodic (free-space) Ewald sum for the electrostatic potential,
%using the Spectral Ewald method.
%
%   out = SE0P_Laplace_fourier_space(x, f, opt)
%   out = SE0P_Laplace_fourier_space(x, f, opt, pre_kernel)
%
%   Computes the electrostatic potential from N point charges.
%
%   Input parameters:
%       :param x: source locations (N×3)
%       :param f: source charges (N×1)
%       :param opt: option struct, see "help SE0P_base_fourier_space"
%       :param pre_kernel: optional, output from SE0P_Laplace_precompute_kernel_fft
%
%   :returns: **out** -- output struct, see "help SE0P_base_fourier_space"

if nargin < 4 || isempty(pre_kernel)
  pre_kernel = SE0P_Laplace_precompute_kernel_fft(opt);
end

opt.kernel = 'laplace';
out = SE0P_base_fourier_space(x, f, opt, pre_kernel);

end
