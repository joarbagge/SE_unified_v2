function out = SE2P_Stokeslet_fourier_space(x, f, opt)
%SE2P_Stokeslet_fourier_space  Compute Fourier-space part of the
%2-periodic Ewald sum for the stokeslet potential, using the
%Spectral Ewald method.
%
%   out = SE2P_Stokeslet_fourier_space(x, f, opt)
%
%   Computes the Stokesian velocity field from N point forces.
%   Periodicity is assumed in the first two directions (x and y),
%   while the third direction (z) is assumed free.
%
%   Input parameters:
%       :param x: source locations (N×3)
%       :param f: source forces (N×3)
%       :param opt: option struct, see "help SE2P_base_fourier_space"
%
%   :returns: **out** -- output struct, see "help SE2P_base_fourier_space"

opt.kernel = 'stokeslet';
if ~isfield(opt,'source_quantity'), opt.source_quantity = sum(f(:).^2); end
out = SE2P_base_fourier_space(x, f, opt);

end
