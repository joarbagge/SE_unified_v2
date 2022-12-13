function H = SE0P_Laplace_scaling(H, opt, pre_kernel, pre_window)
%SE0P_Laplace_scaling  Fourier-space scaling for 0-periodic (free-space) electrostatics
%
%   H = SE0P_Laplace_scaling(H, opt, pre_kernel, pre_window)
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
TMP = K1.^2 + K2.^2 + K3.^2;
if ~((opt.return_fouriervalgrad || opt.return_gridvalgrad || opt.return_potgrad) ...
    && opt.use_ik_diff)
  clear K1 K2 K3;
end

% The scaling is given by: Green × Screening × Window^(-2),
% each factor representing the Fourier transform of a function.
% For the Gaussian window, we combine Screening×Window^(-2) into
% a single factor (Combo).

TMP = TMP/(4*opt.xi^2);
if strcmp(opt.window, 'gaussian') % special treatment for Gaussian window
  eta = 2*(w*opt.xi)^2 / opt.window_shape;
  % Screening and Window factors combined
  % 1. Using Fourier transform of untruncated Gaussian
  TMP = exp(-TMP*(1-eta*pw/2)); % Combo
  TMP = pre_kernel.kernel .* TMP; % multiply Combo by Green to get Scaling
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
  TMP = pre_kernel.kernel .* Window_m2 .* TMP; % Scaling
  clear Window_m2;
end
H{1} = TMP .* H{1};

% ik differentiation to compute gradient
if (opt.return_fouriervalgrad || opt.return_gridvalgrad || opt.return_potgrad) ...
    && opt.use_ik_diff
  Hi = 1i * H{1};
  H{2} = Hi .* K1;
  H{3} = Hi .* K2;
  H{4} = Hi .* K3;
end
if ~(opt.return_fourierval || opt.return_gridval || opt.return_pot ...
    || (opt.return_potgrad && ~opt.use_ik_diff))
  H{1} = [];
end

% ------------------------------------------------------------------------------
function F = kaiser_exact_ft(k2, b2, w, scaling)

t = sqrt(b2 - k2*w^2);
F = 2*w*sinh(t)./t * scaling;
