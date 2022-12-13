function out = SE3P_Stresslet_fourier_space(x, q, n, opt)
%SE3P_Stresslet_fourier_space  Compute Fourier-space part of the
%3-periodic Ewald sum for the stresslet potential, using the
%Spectral Ewald method.
%
%   out = SE3P_Stresslet_fourier_space(x, q, n, opt)
%
%   Computes the stresslet velocity field from N point sources.
%
%   Input parameters:
%       :param x: source locations (N×3)
%       :param q: source strengths (N×3)
%       :param n: source normals (N×3)
%       :param opt: option struct, see "help SE3P_base_fourier_space"
%
%   :returns: **out** -- output struct, see "help SE3P_base_fourier_space"

N = size(x, 1);
f = zeros(N, 3, 3);
for j=1:3
  for l=1:3
    f(:,j,l) = q(:,j) .* n(:,l);
  end
end

opt.kernel = 'stresslet';
out = SE3P_base_fourier_space(x, f, opt);

end
