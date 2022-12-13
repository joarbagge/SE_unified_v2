function [H, Hl, H0] = SE1P_Rotlet_scaling(H, Hl, H0, opt, pre_fft)
%SE1P_Rotlet_scaling  Fourier-space scaling for 1-periodic rotlet
%
%   [H, Hl, H0] = SE1P_Rotlet_scaling(H, Hl, H0, opt, pre_fft)
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
K = cell(3,1);
k1 = SE_aft_k_vectors(opt.grid(1), opt.box(1), 1);
k2 = SE_aft_k_vectors(opt.grid2, opt.box2, opt.actual_upsampling_global(1)); % free direction
k3 = SE_aft_k_vectors(opt.grid3, opt.box3, opt.actual_upsampling_global(2)); % free direction
[K{1}, K{2}, K{3}] = ndgrid(k1, k2, k3);
KK = K{1}.^2 + K{2}.^2 + K{3}.^2;

% The scaling is given by: Green × Screening × Window^(-2),
% each factor representing the Fourier transform of a function.
% For the Gaussian window, we combine Screening×Window^(-2) into
% a single factor (Combo).

% For the rotlet, we first use the Laplace Green's function,
% and then use a relation between it and the rotlet.

% Scaling for the global Fourier domain
Laplace = 4*pi./KK;
C = KK/(4*opt.xi^2);
if strcmp(opt.window, 'gaussian') % special treatment for Gaussian window
  eta = 2*(w*opt.xi)^2 / opt.window_shape;
  % Screening and Window factors combined
  Combo = exp(-C*(1-eta*pw/2));
  Scaling = 1i * Laplace .* Combo;
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
  Scaling = 1i * Laplace .* Screening .* Window_m2;
end
% Now the element corresponding to the zero mode is Inf.
% This is fine since we overwrite it later on.
if isempty(H0{1})
  H0_cache = cell(size(H));
  for j=1:numel(H)
    H0_cache{j} = squeeze(H{j}(1,1,:));
  end
end
% Use relation between rotlet and Laplace Green's function
KcrossH = cell(3,1);
KcrossH{1} = K{2}.*H{3} - K{3}.*H{2};
KcrossH{2} = K{3}.*H{1} - K{1}.*H{3};
KcrossH{3} = K{1}.*H{2} - K{2}.*H{1};
for j=1:3
  H{j} = Scaling .* KcrossH{j};
end

% Scaling for the local pad (same structure as above)
if numel(opt.local_modes1) > 0
  n1 = opt.local_modes1;
  k2 = SE_aft_k_vectors(opt.grid2, opt.box2, opt.actual_upsampling_local(1));
  k3 = SE_aft_k_vectors(opt.grid3, opt.box3, opt.actual_upsampling_local(2));
  [K{1}, K{2}, K{3}] = ndgrid(k1(n1), k2, k3);
  KK = K{1}.^2 + K{2}.^2 + K{3}.^2;
  Laplace = 4*pi./KK;
  C = KK/(4*opt.xi^2);
  if strcmp(opt.window, 'gaussian')
    Combo = exp(-C*(1-eta*pw/2));
    Scaling_local = 1i * Laplace .* Combo;
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
    Scaling_local = 1i * Laplace .* Screening .* Window_m2;
  end
  % Use relation between rotlet and Laplace Green's function
  KcrossH = cell(3,1);
  KcrossH{1} = K{2}.*Hl{3} - K{3}.*Hl{2};
  KcrossH{2} = K{3}.*Hl{1} - K{1}.*Hl{3};
  KcrossH{3} = K{1}.*Hl{2} - K{2}.*Hl{1};
  for j=1:3
    Hl{j} = Scaling_local .* KcrossH{j};
  end
end

% Scaling for the zero mode [k1==0] (Vico method)
K{1} = 0;
k2 = SE_aft_k_vectors(opt.grid2, opt.box2, opt.actual_upsampling_zero(1));
k3 = SE_aft_k_vectors(opt.grid3, opt.box3, opt.actual_upsampling_zero(2));
[K{2}, K{3}] = ndgrid(k2, k3);
KK = K{2}.^2 + K{3}.^2;
C = KK/(4*opt.xi^2);
Kn = sqrt(KK);
R = opt.greens_truncation_R;
% Modified Green's function (based on 4*pi/k^2)
if ~isfield(opt, 'vico_variant') || opt.vico_variant == 1
  %lH = R; % optimal value for Vico formula convergence
  %Green = 4*pi*((1-besselj(0,R*Kn))./KK - R*log(R/lH)*besselj(1,R*Kn)./Kn);
  %Green(1,1) = pi*R^2 * (1 - 2*log(R/lH)); % finite limit at k2==k3==0
  % Use simplified formula, better speed
  Green = 4*pi*(1-besselj(0,R*Kn))./KK;
  Green(1,1) = pi*R^2; % finite limit at k2==k3==0
  Kfactor = -1i * Kn;
elseif opt.vico_variant == 2
  Green = 4*pi*( ...
                 R*besselj(0,R*Kn) ...
                 + pi*R/2*besselj(1,R*Kn).*StruveH0(R*Kn) ...
                 - pi*R/2*besselj(0,R*Kn).*StruveH1(R*Kn) ...
               );
  Green(1,1) = 4*pi*R; % finite limit at k2==k3==0
  Kfactor = -1i;
end
if strcmp(opt.window, 'gaussian')
  Combo = exp(-C*(1-eta*pw/2));
  Scaling_zero = Kfactor .* Green .* Combo;
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
  Scaling_zero = Kfactor .* Green .* Screening .* Window_m2;
end
%% Plot stuff... This should be radial now, right?
%if opt.vico_variant == 2
%  Kfactor = Kfactor*ones(size(K{2}));
%end
%N = numel(k2) / 2;
%clf
%subplot(2,2,1)
%loglog(k2(1:N), abs(Kfactor(1:N,1)), '.-')
%hold on
%loglog(k2(1:N), k2(1:N), 'k--')
%loglog(k2(1:N), ones(1,N), 'k--')
%subplot(2,2,2)
%loglog(k2(1:N), abs(Green(1:N,1)), '.-')
%hold on
%if opt.vico_variant == 1
%  loglog(k2(1:N), 20./k2(1:N).^2, 'k--')
%  loglog(k2(1:N), 20./k2(1:N).^1.5, 'k--')
%elseif opt.vico_variant == 2
%  loglog(k2(1:N), 15./k2(1:N), 'k--')
%end
%subplot(2,2,3)
%loglog(k2(1:N), abs(Kfactor(1:N,1) .* Green(1:N,1)), '.-')
%hold on
%if opt.vico_variant == 1
%  loglog(k2(1:N), 20./k2(1:N).^1, 'k--')
%  loglog(k2(1:N), 20./k2(1:N).^0.5, 'k--')
%elseif opt.vico_variant == 2
%  loglog(k2(1:N), 15./k2(1:N), 'k--')
%end
if ~isempty(H0{1})
  H0_cache = H0;
end
% Use simplified relation between rotlet and Laplace Green's function
K2n = K{2} ./ Kn;
K3n = K{3} ./ Kn;
K2n(1,1) = 0;
K3n(1,1) = 0;
%HH = H0_cache;
tmp = H0_cache{1};
H0_cache{1} = Scaling_zero .* (H0_cache{2}.*K3n - H0_cache{3}.*K2n);
H0_cache{2} = -Scaling_zero .* tmp.*K3n;
H0_cache{3} = Scaling_zero .* tmp.*K2n;
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
