function out = SE1P_Stokeslet_fourier_space(x, f, opt)
%SE1P_Stokeslet_fourier_space  Compute Fourier-space part of the
%1-periodic Ewald sum for the stokeslet potential, using the
%Spectral Ewald method.
%
%   out = SE1P_Stokeslet_fourier_space(x, f, opt)
%
%   Computes the Stokesian velocity field from N point forces.
%   Periodicity is assumed in the first direction (x), while the
%   two remaining directions (y and z) are assumed free.
%
%   Input parameters:
%       :param x: source locations (N×3)
%       :param f: source forces (N×3)
%       :param opt: option struct, see "help SE1P_base_fourier_space"
%
%   The stokeslet has an additional option controlling the choice
%   of constant in the Fourier space zero mode:
%       :param opt.stokeslet_k0_constant: the value of c in the
%           2D biharmonic Green's function -r^2 (log(r) - c).
%           (default: 0)
%
%   :returns: **out** -- output struct, see "help SE1P_base_fourier_space"

opt.kernel = 'stokeslet';
if ~isfield(opt,'source_quantity'), opt.source_quantity = sum(f(:).^2); end
out = SE1P_base_fourier_space(x, f, opt);

end
