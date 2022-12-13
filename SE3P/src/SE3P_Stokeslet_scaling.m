function H = SE3P_Stokeslet_scaling(H, opt, pre_fft)
%SE3P_Stokeslet_scaling  Fourier-space scaling for 3-periodic stokeslet
%
%   H = SE3P_Stokeslet_scaling(H, opt, pre_fft)
%
%   Input parameters:
%       :param H: input data in Fourier space (cell array)
%       :param opt: structure with Ewald options
%       :param pre_fft: optional precomputation structure:
%       :param pre_fft.window: W^(-pw) where W is the Fourier transform of the window function
%
%   :returns: **H** -- output data in Fourier space (cell array)

if nargin < 3
  pre_fft = [];
end
pw = opt.window_scaling_power; % should be 2 or 1 typically
w = opt.window_halfwidth;

% Compute k-vectors
[k1, k2, k3] = SE_k_vectors(opt.grid, opt.box);
[K1, K2, K3] = ndgrid(k1, k2, k3);
KK = K1.^2 + K2.^2 + K3.^2;

% The scaling is given by: Green × Screening × Window^(-2),
% each factor representing the Fourier transform of a function.
% For the Gaussian window, we combine Screening×Window^(-2) into
% a single factor (Combo).

% For the stokeslet, we first use the biharmonic Green's function,
% and then use a relation between it and the stokeslet.

Biharmonic = 8*pi./(KK.*KK);
C = KK/(4*opt.xi^2);
if strcmp(opt.window, 'gaussian')
  eta = 2*(w*opt.xi)^2 / opt.window_shape;
  % Screening and Window factors combined
  Combo = (1+C) .* exp(-C*(1-eta*pw/2));
  Scaling = Biharmonic .* Combo;
  % FIXME: LEGACY code used the MEX function below (needs complex H):
  %[H1 H2 H3] = se3p_fast_k_scaling(H1,H2,H3,xi,opt.box,eta);
else
  Screening = (1+C) .* exp(-C);
  if ~isempty(pre_fft)
    Window_m2 = pre_fft.window;
  elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
    % Compute Window_m2 on the fly
    b2 = opt.window_shape^2;
    f1 = kaiser_exact_ft(k1.^2, b2, w, opt.kaiser_scaling);
    f2 = kaiser_exact_ft(k2.^2, b2, w, opt.kaiser_scaling);
    f3 = kaiser_exact_ft(k3.^2, b2, w, opt.kaiser_scaling);
    [F1, F2, F3] = ndgrid(f1.^pw, f2.^pw, f3.^pw);
    F = F1.*F2.*F3; % tensor product of spatial directions
    Window_m2 = 1./F;
  else
    error('Precomputed Fourier transform necessary for the %s window', opt.window);
  end
  Scaling = Biharmonic .* Screening .* Window_m2;
end
Scaling(1,1,1) = 0; % zero mode should be zero

% Use relation between stokeslet and biharmonic Green's function
KdotH = K1.*H{1} + K2.*H{2} + K3.*H{3};
H{1} = Scaling .* (KK.*H{1} - KdotH.*K1);
H{2} = Scaling .* (KK.*H{2} - KdotH.*K2);
H{3} = Scaling .* (KK.*H{3} - KdotH.*K3);

% ------------------------------------------------------------------------------
function F = kaiser_exact_ft(k2, b2, w, scaling)

t = sqrt(b2 - k2*w^2);
F = 2*w*sinh(t)./t * scaling;
