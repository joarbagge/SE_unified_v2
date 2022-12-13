function out = SE3P_Stokeslet_fourier_space(x, f, opt)
%SE3P_Stokeslet_fourier_space  Compute Fourier-space part of the
%3-periodic Ewald sum for the stokeslet potential, using the
%Spectral Ewald method.
%
%   out = SE3P_Stokeslet_fourier_space(x, f, opt)
%
%   Computes the Stokesian velocity field from N point forces.
%
%   Input parameters:
%       :param x: source locations (N×3)
%       :param f: source forces (N×3)
%       :param opt: option struct, see "help SE3P_base_fourier_space"
%
%   :returns: **out** -- output struct, see "help SE3P_base_fourier_space"

opt.kernel = 'stokeslet';
out = SE3P_base_fourier_space(x, f, opt);

end
