function [H, Hl, H0] = SE1P_Laplace_scaling(H, Hl, H0, opt, pre_fft)
%SE1P_Laplace_scaling  Fourier-space scaling for 1-periodic electrostatics
%
%   [H, Hl, H0] = SE1P_Laplace_scaling(H, Hl, H0, opt, pre_fft)
%
%   Input parameters:
%       :param [H, Hl, H0]: input data in Fourier space (cell arrays)
%       :param opt: structure with Ewald options
%       :param pre_fft: optional precomputation structure:
%       :param pre_fft.window: W^(-pw) where W is the Fourier transform of the window function
%       :param pre_fft.window_local: W^(-pw) on the local pad
%       :param pre_fft.window_zero: W^(-pw) for the zero mode
%
%   :returns: **[H, Hl, H0]** -- output data in Fourier space (cell arrays)

if nargin < 5
  pre_fft = [];
end
pw = opt.window_scaling_power; % should be 2 or 1 typically
w = opt.window_halfwidth;

% Compute k-vectors
k1 = SE_aft_k_vectors(opt.grid(1), opt.box(1), 1);
k2 = SE_aft_k_vectors(opt.grid2, opt.box2, opt.actual_upsampling_global(1)); % free direction
k3 = SE_aft_k_vectors(opt.grid3, opt.box3, opt.actual_upsampling_global(2)); % free direction
[K1, K2, K3] = ndgrid(k1, k2, k3);
KK = K1.^2 + K2.^2 + K3.^2;

% The scaling is given by: Green × Screening × Window^(-2),
% each factor representing the Fourier transform of a function.
% For the Gaussian window, we combine Screening×Window^(-2) into
% a single factor (Combo).

% For Laplace, Green = 4*pi/KK.

% Scaling for the global Fourier domain
C = KK/(4*opt.xi^2);
if strcmp(opt.window, 'gaussian') % special treatment for Gaussian window
  eta = 2*(w*opt.xi)^2 / opt.window_shape;
  % Screening and Window factors combined
  % 1. Using Fourier transform of untruncated Gaussian
  Combo = exp(-C*(1-eta*pw/2));
  Scaling = 4*pi * Combo ./ KK; % multiply by Green
else
  Screening = exp(-C);
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
  Scaling = 4*pi * Window_m2 .* Screening ./ KK;
end
% Now the element corresponding to the zero mode is Inf.
% This is fine since we overwrite it later on.
if isempty(H0{1})
  H0_cache = cell(size(H));
  for j=1:numel(H)
    H0_cache{j} = squeeze(H{j}(1,:,:));
  end
end
H{1} = Scaling .* H{1};
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

% Scaling for the local pad (same structure as above)
if numel(opt.local_modes1) > 0
  n1 = opt.local_modes1;
  k2 = SE_aft_k_vectors(opt.grid2, opt.box2, opt.actual_upsampling_local(1));
  k3 = SE_aft_k_vectors(opt.grid3, opt.box3, opt.actual_upsampling_local(2));
  [K1, K2, K3] = ndgrid(k1(n1), k2, k3);
  KK = K1.^2 + K2.^2 + K3.^2;
  C = KK/(4*opt.xi^2);
  if strcmp(opt.window, 'gaussian')
    Combo = exp(-C*(1-eta*pw/2));
    Scaling_local = 4*pi * Combo ./ KK;
  else
    Screening = exp(-C);
    if ~isempty(pre_fft)
      Window_m2 = pre_fft.window_local;
    elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
      f2 = kaiser_exact_ft(k2.^2, b2, w, opt.kaiser_scaling);
      f3 = kaiser_exact_ft(k3.^2, b2, w, opt.kaiser_scaling);
      [F1, F2, F3] = ndgrid(f1(n1).^pw, f2.^pw, f3.^pw);
      F = F1.*F2.*F3;
      Window_m2 = 1./F;
    else
      error('Precomputed Fourier transform necessary for the %s window', opt.window);
    end
    Scaling_local = 4*pi * Window_m2 .* Screening ./ KK;
  end
  Hl{1} = Scaling_local .* Hl{1};
  % ik differentiation to compute gradient
  if (opt.return_fouriervalgrad || opt.return_gridvalgrad || opt.return_potgrad) ...
      && opt.use_ik_diff
    Hi = 1i * Hl{1};
    Hl{2} = Hi .* K1;
    Hl{3} = Hi .* K2;
    Hl{4} = Hi .* K3;
  end
  if ~(opt.return_fourierval || opt.return_gridval || opt.return_pot ...
      || (opt.return_potgrad && ~opt.use_ik_diff))
    Hl{1} = [];
  end
end

% Scaling for the zero mode [k1==0] (Vico method)
k2 = SE_aft_k_vectors(opt.grid2, opt.box2, opt.actual_upsampling_zero(1));
k3 = SE_aft_k_vectors(opt.grid3, opt.box3, opt.actual_upsampling_zero(2));
[K2, K3] = ndgrid(k2, k3);
KK = k1(1).^2 + K2.^2 + K3.^2;
C = KK/(4*opt.xi^2);
K = sqrt(KK);
R = opt.greens_truncation_R;
% Modified Green's function
Green = 4*pi*((1-besselj(0,R*K))./KK - R*log(R)*besselj(1,R*K)./K);
Green(1,1) = 2*pi*R^2 * (1/2 - log(R)); % finite limit at k2==k3==0
if strcmp(opt.window, 'gaussian')
  Combo = exp(-C*(1-eta*pw/2));
  Scaling_zero = Green .* Combo;
else
  Screening = exp(-C);
  if ~isempty(pre_fft)
    Window_m2 = pre_fft.window_zero;
  elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
    f2 = kaiser_exact_ft(k2.^2, b2, w, opt.kaiser_scaling);
    f3 = kaiser_exact_ft(k3.^2, b2, w, opt.kaiser_scaling);
    F1 = f1(1).^pw; [F2, F3] = ndgrid(f2.^pw, f3.^pw);
    F = F1.*F2.*F3;
    Window_m2 = 1./F;
  else
    error('Precomputed Fourier transform necessary for the %s window', opt.window);
  end
  Scaling_zero = Green .* Window_m2 .* Screening;
end
if ~isempty(H0{1})
  H0_cache = H0;
end
H0_cache{1} = Scaling_zero .* H0_cache{1};
% ik differentiation to compute gradient
if (opt.return_fouriervalgrad || opt.return_gridvalgrad || opt.return_potgrad) ...
    && opt.use_ik_diff
  Hi = 1i * H0_cache{1};
  H0_cache{2} = Hi .* k1(1);
  H0_cache{3} = Hi .* K2;
  H0_cache{4} = Hi .* K3;
end
if ~(opt.return_fourierval || opt.return_gridval || opt.return_pot ...
    || (opt.return_potgrad && ~opt.use_ik_diff))
  H0_cache{1} = [];
end
if isempty(H0{1})
  % Store in H instead
  for j=1:numel(H)
    H{j}(1,:,:) = H0_cache{j};
  end
else
  H0 = H0_cache;
end

% ------------------------------------------------------------------------------
function F = kaiser_exact_ft(k2, b2, w, scaling)

t = sqrt(b2 - k2*w^2);
F = 2*w*sinh(t)./t * scaling;
