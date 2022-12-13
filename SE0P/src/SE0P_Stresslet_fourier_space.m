function out = SE0P_Stresslet_fourier_space(x, q, n, opt, pre_kernel)
%SE0P_Stresslet_fourier_space  Compute Fourier-space part of the
%0-periodic (free-space) Ewald sum for the stresslet potential,
%using the Spectral Ewald method.
%
%   out = SE0P_Stresslet_fourier_space(x, q, n, opt)
%   out = SE0P_Stresslet_fourier_space(x, q, n, opt, pre_kernel)
%
%   Computes the stresslet velocity field from N point sources.
%
%   Input parameters:
%       :param x: source locations (N×3)
%       :param q: source strengths (N×3)
%       :param n: source normals (N×3)
%       :param opt: option struct, see "help SE0P_base_fourier_space"
%       :param pre_kernel: optional, output from SE0P_Stresslet_precompute_kernel_fft
%
%   :returns: **out** -- output struct, see "help SE0P_base_fourier_space"

if nargin < 5 || isempty(pre_kernel)
  pre_kernel = SE0P_Stresslet_precompute_kernel_fft(opt);
end

N = size(x, 1);
f = zeros(N, 3, 3);
for j=1:3
  for l=1:3
    f(:,j,l) = q(:,j) .* n(:,l);
  end
end

opt.kernel = 'stresslet';
out = SE0P_base_fourier_space(x, f, opt, pre_kernel);

end
