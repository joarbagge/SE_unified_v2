function H = SE0P_Rotlet_scaling(H, opt, pre_kernel, pre_window)
%SE0P_Rotlet_scaling  Fourier-space scaling for 0-periodic (free-space) rotlet
%
%   H = SE0P_Rotlet_scaling(H, opt, pre_kernel, pre_window)
%
%   Input parameters:
%       :param H: input data in Fourier space (cell array)
%       :param opt: structure with Ewald options
%       :param pre_kernel: mandatory precomputation structure:
%       :param pre_kernel.kernel: Fourier transform of the harmonic Green's function
%       :param pre_window: optional precomputation structure:
%       :param pre_window.window: W^(-pw) where W is the Fourier transform of the window function
%
%   :returns: **H** -- output data in Fourier space (cell array)

if nargin < 4
  pre_window = [];
end
pw = opt.window_scaling_power; % should be 2 or 1 typically
w = opt.window_halfwidth;

% Compute k-vectors
[k1, k2, k3] = SE_k_vectors(opt.padded_grid, opt.padded_box, 'shifted');
[K1, K2, K3] = ndgrid(k1, k2, k3);
KK = K1.^2 + K2.^2 + K3.^2;

% The scaling is given by: Green × Screening × Window^(-2),
% each factor representing the Fourier transform of a function.
% For the Gaussian window, we combine Screening×Window^(-2) into
% a single factor (Combo).

% For the rotlet, we first use the harmonic Green's function,
% and then use a relation between it and the rotlet.

Green = pre_kernel.kernel;
if ~isfield(opt, 'vico_variant') || opt.vico_variant == 1
  Kfactor = 1i;
elseif opt.vico_variant == 2
  Kfactor = -0.5i * KK;
end
TMP = KK/(4*opt.xi^2);
if strcmp(opt.window, 'gaussian') % special treatment for Gaussian window
  eta = 2*(w*opt.xi)^2 / opt.window_shape;
  % Screening and Window factors combined
  TMP = exp(-TMP*(1-eta*pw/2)); % Combo
  TMP = Kfactor .* Green .* TMP; % multiply Combo by Green to get Scaling
else
  TMP = exp(-TMP); % Screening
  if ~isempty(pre_window)
    Window_m2 = pre_window.window;
  elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
    % Compute Window_m2 on the fly
    b2 = opt.window_shape^2;
    f1 = kaiser_exact_ft(k1.^2, b2, w, opt.kaiser_scaling);
    f2 = kaiser_exact_ft(k2.^2, b2, w, opt.kaiser_scaling);
    f3 = kaiser_exact_ft(k3.^2, b2, w, opt.kaiser_scaling);
    [F1, F2, F3] = ndgrid(f1.^pw, f2.^pw, f3.^pw);
    Window_m2 = 1./(F1.*F2.*F3); % tensor product of spatial directions
    clear F1 F2 F3;
  else
    error('Precomputed Fourier transform necessary for the %s window', opt.window);
  end
  TMP = Kfactor .* Green .* TMP .* Window_m2; % Scaling
  clear Window_m2;
end
clear Green Kfactor;
% Use relation between rotlet and harmonic Green's function
H1 = H{1}; H2 = H{2};
H{1} = TMP .* (K2.*H{3} - K3.*H{2});
H{2} = TMP .* (K3.*H1   - K1.*H{3});
H{3} = TMP .* (K1.*H2   - K2.*H1  );

% ------------------------------------------------------------------------------
function F = kaiser_exact_ft(k2, b2, w, scaling)

t = sqrt(b2 - k2*w^2);
F = 2*w*sinh(t)./t * scaling;
