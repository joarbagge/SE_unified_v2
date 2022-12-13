function out = SE3P_Laplace_fourier_space(x, f, opt)
%SE3P_Laplace_fourier_space  Compute Fourier-space part of the
%3-periodic Ewald sum for the electrostatic potential, using the
%Spectral Ewald method.
%
%   out = SE3P_Laplace_fourier_space(x, f, opt)
%
%   Computes the electrostatic potential from N point charges.
%
%   Input parameters:
%       :param x: source locations (N×3)
%       :param f: source charges (N×1)
%       :param opt: option struct, see "help SE3P_base_fourier_space"
%
%   :returns: **out** -- output struct, see "help SE3P_base_fourier_space"

opt.kernel = 'laplace';
out = SE3P_base_fourier_space(x, f, opt);

end
