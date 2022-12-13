function [H, Hl, H0] = SE2P_Stokeslet_scaling(H, Hl, H0, opt, pre_fft)
%SE2P_Stokeslet_scaling  Fourier-space scaling for 2-periodic stokeslet
%
%   [H, Hl, H0] = SE2P_Stokeslet_scaling(H, Hl, H0, opt, pre_fft)
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
k2 = SE_aft_k_vectors(opt.grid(2), opt.box(2), 1);
k3 = SE_aft_k_vectors(opt.grid3, opt.box3, opt.actual_upsampling_global); % free direction
[K1, K2, K3] = ndgrid(k1, k2, k3);
KK = K1.^2 + K2.^2 + K3.^2;

% The scaling is given by: Green × Screening × Window^(-2),
% each factor representing the Fourier transform of a function.
% For the Gaussian window, we combine Screening×Window^(-2) into
% a single factor (Combo).

% For the stokeslet, we first use the biharmonic Green's function,
% and then use a relation between it and the stokeslet.

% Scaling for the global Fourier domain
Biharmonic = 8*pi./(KK.*KK);
C = KK/(4*opt.xi^2);
if strcmp(opt.window, 'gaussian') % special treatment for Gaussian window
  eta = 2*(w*opt.xi)^2 / opt.window_shape;
  % Screening and Window factors combined
  Combo = (1+C) .* exp(-C*(1-eta*pw/2));
  Scaling = Biharmonic .* Combo;
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
% Now the element corresponding to the zero mode is Inf.
% This is fine since we overwrite it later on.
if isempty(H0{1})
  H0_cache = cell(size(H));
  for j=1:numel(H)
    H0_cache{j} = squeeze(H{j}(1,1,:));
  end
end
% Use relation between stokeslet and biharmonic Green's function
KdotH = K1.*H{1} + K2.*H{2} + K3.*H{3};
H{1} = Scaling .* (KK.*H{1} - KdotH.*K1);
H{2} = Scaling .* (KK.*H{2} - KdotH.*K2);
H{3} = Scaling .* (KK.*H{3} - KdotH.*K3);

% Scaling for the local pad (same structure as above)
if numel(opt.local_modes1) > 0 || numel(opt.local_modes2) > 0
  n1 = opt.local_modes1;
  n2 = opt.local_modes2;
  k3 = SE_aft_k_vectors(opt.grid3, opt.box3, opt.actual_upsampling_local);
  [K1, K2, K3] = ndgrid(k1(n1), k2(n2), k3);
  KK = K1.^2 + K2.^2 + K3.^2;
  Biharmonic = 8*pi./(KK.*KK);
  C = KK/(4*opt.xi^2);
  if strcmp(opt.window, 'gaussian')
    Combo = (1+C) .* exp(-C*(1-eta*pw/2));
    Scaling_local = Biharmonic .* Combo;
  else
    Screening = (1+C) .* exp(-C);
    if ~isempty(pre_fft)
      Window_m2 = pre_fft.window_local;
    elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
      f3 = kaiser_exact_ft(k3.^2, b2, w, opt.kaiser_scaling);
      [F1, F2, F3] = ndgrid(f1(n1).^pw, f2(n2).^pw, f3.^pw);
      F = F1.*F2.*F3;
      Window_m2 = 1./F;
    else
      error('Precomputed Fourier transform necessary for the %s window', opt.window);
    end
    Scaling_local = Biharmonic .* Screening .* Window_m2;
  end
  % Use relation between stokeslet and biharmonic Green's function
  KdotH = K1.*Hl{1} + K2.*Hl{2} + K3.*Hl{3};
  Hl{1} = Scaling_local .* (KK.*Hl{1} - KdotH.*K1);
  Hl{2} = Scaling_local .* (KK.*Hl{2} - KdotH.*K2);
  Hl{3} = Scaling_local .* (KK.*Hl{3} - KdotH.*K3);
end

% Scaling for the zero mode [k1==k2==0] (Vico method)
k3 = SE_aft_k_vectors(opt.grid3, opt.box3, opt.actual_upsampling_zero)';
KK = k1(1).^2 + k2(1).^2 + k3.^2;
C = KK/(4*opt.xi^2);
Kn = sqrt(KK);
R = opt.greens_truncation_R;
% Modified Green's function
if ~isfield(opt, 'vico_variant') || opt.vico_variant == 3
  Green = -8*pi * (1 - cos(R*Kn) - R*Kn.*sin(R*Kn)) ./ KK;
  Green(1) = 4*pi*R^2; % finite limit at k3==0
  Kfactor = -1;
elseif opt.vico_variant == 2
  Green = -4*pi*( ...
                  -2*R*Kn.*sin(R*Kn) ...
                  + (R^2*KK-2).*cos(R*Kn) ...
                  + 2 ...
                ) ./ (KK.*k3);
  Green(1) = 0; % finite limit at k3==0
  Kfactor = -k3;
elseif opt.vico_variant == 1
  %Green = -4*pi*( ...
  %                R*Kn.*(R^2*KK-6).*sin(R*Kn) ...
  %                + 3*(R^2*KK-2).*cos(R*Kn) ...
  %                + 6 ...
  %              ) ./ (3*KK.*KK);
  %Green = (4/3)*pi*R*Kn.*(R^2-6./KK).*sin(R*Kn) ...
  %        + (4/3)*pi*3*(R^2-2./KK).*cos(R*Kn) ...
  %        + (4/3)*pi*6./KK;
  %Green = (4/3)*pi*R*Kn.*(R^2-6./KK).*sin(R*Kn);
  %Green = (4/3)*pi*3*(R^2-2./KK).*cos(R*Kn) + (4/3)*pi*6./KK;
  Green = 4*pi*(R^2*cos(R*Kn) + 2*(1-cos(R*Kn))./KK);
  %Green(1) = -pi*R^4/3; % finite limit at k3==0
  Green(1) = 0; % finite limit at k3==0
  Kfactor = 1;
end
if strcmp(opt.window, 'gaussian')
  Combo = (1+C) .* exp(-C*(1-eta*pw/2));
  Scaling_zero = Kfactor .* Green .* Combo;
else
  Screening = (1+C) .* exp(-C);
  if ~isempty(pre_fft)
    Window_m2 = pre_fft.window_zero;
  elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
    f3 = kaiser_exact_ft(k3.^2, b2, w, opt.kaiser_scaling);
    F1 = f1(1).^pw; F2 = f2(1).^pw; F3 = f3.^pw;
    F = F1*F2*F3;
    Window_m2 = 1./F;
  else
    error('Precomputed Fourier transform necessary for the %s window', opt.window);
  end
  Scaling_zero = Kfactor .* Green .* Screening .* Window_m2;
end
if ~isempty(H0{1})
  H0_cache = H0;
end
H0_cache{1} = Scaling_zero .* H0_cache{1};
H0_cache{2} = Scaling_zero .* H0_cache{2};
H0_cache{3} = 0*H0_cache{3};
if isempty(H0{1})
  % Store in H instead
  for j=1:numel(H)
    H{j}(1,1,:) = H0_cache{j};
  end
else
  H0 = H0_cache;
end

% ------------------------------------------------------------------------------
function F = kaiser_exact_ft(k2, b2, w, scaling)

t = sqrt(b2 - k2*w^2);
F = 2*w*sinh(t)./t * scaling;
