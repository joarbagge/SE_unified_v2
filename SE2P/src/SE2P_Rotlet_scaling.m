function [H, Hl, H0] = SE2P_Rotlet_scaling(H, Hl, H0, opt, pre_fft)
%SE2P_Rotlet_scaling  Fourier-space scaling for 2-periodic rotlet
%
%   [H, Hl, H0] = SE2P_Rotlet_scaling(H, Hl, H0, opt, pre_fft)
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
k2 = SE_aft_k_vectors(opt.grid(2), opt.box(2), 1);
k3 = SE_aft_k_vectors(opt.grid3, opt.box3, opt.actual_upsampling_global); % free direction
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
if numel(opt.local_modes1) > 0 || numel(opt.local_modes2) > 0
  n1 = opt.local_modes1;
  n2 = opt.local_modes2;
  k3 = SE_aft_k_vectors(opt.grid3, opt.box3, opt.actual_upsampling_local);
  [K{1}, K{2}, K{3}] = ndgrid(k1(n1), k2(n2), k3);
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
      f3 = kaiser_exact_ft(k3.^2, b2, w, opt.kaiser_scaling);
      [F1, F2, F3] = ndgrid(f1(n1).^pw, f2(n2).^pw, f3.^pw);
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

% Scaling for the zero mode [k1==k2==0] (Vico method)
K{1} = 0;
K{2} = 0;
k3 = SE_aft_k_vectors(opt.grid3, opt.box3, opt.actual_upsampling_zero)';
K{3} = k3;
KK = K{1}.^2 + K{2}.^2 + K{3}.^2;
C = KK/(4*opt.xi^2);
Kn = sqrt(KK);
R = opt.greens_truncation_R;
% Modified Green's function (based on 4*pi/k3)
if ~isfield(opt, 'vico_variant') || opt.vico_variant == 2
  Green = 4*pi * (1 - cos(R*Kn)) ./ k3;
  Green(1) = 0; % finite limit at k3==0
  Kfactor = -1i;
elseif opt.vico_variant == 1
  Green = 4*pi * (1 - cos(R*Kn) - R*Kn.*sin(R*Kn)) ./ KK;
  Green(1) = -2*pi*R^2; % finite limit at k3==0
  Kfactor = -1i*k3;
end
if strcmp(opt.window, 'gaussian')
  Combo = exp(-C*(1-eta*pw/2));
  Scaling_zero = Kfactor .* Green .* Combo;
else
  Screening = exp(-C);
  if ~isempty(pre_fft)
    Window_m2 = pre_fft.window_zero;
  elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
    f3 = kaiser_exact_ft(k3.^2, b2, w, opt.kaiser_scaling);
    F1 = f1(1).^pw; F2 = f2(1).^pw; F3 = f3.^pw;
    F = F1.*F2.*F3;
    Window_m2 = 1./F;
  else
    error('Precomputed Fourier transform necessary for the %s window', opt.window);
  end
  Scaling_zero = Kfactor .* Green .* Screening .* Window_m2;
end
if ~isempty(H0{1})
  H0_cache = H0;
end
% Use simplified relation between rotlet and Laplace Green's function
%HH = H0_cache;
tmp = H0_cache{1};
H0_cache{1} = Scaling_zero .* H0_cache{2};
H0_cache{2} = -Scaling_zero .* tmp;
H0_cache{3} = 0*H0_cache{3};
if isempty(H0{1})
  % Store in H instead
  for j=1:numel(H)
    H{j}(1,1,:) = H0_cache{j};
  end
else
  H0 = H0_cache;
end

% (1), regularizing 4*pi/k^2 (whole harmonic)
%Green_1 = 4*pi * (1 - cos(R*Kn) - R*Kn.*sin(R*Kn)) ./ KK;
%Green_1(1) = -2*pi*R^2; % finite limit at k3==0
%Scaling_zero = -1i .* k3 .* Green_1 .* Screening .* Window_m2;
%H0_1 = cell(3,1);
%H0_1{1} = Scaling_zero .* HH{2};
%H0_1{2} = -Scaling_zero .* HH{1};
%H0_1{3} = 0*HH{3};

% (2), regularizing 4*pi/k
%Green_2 = 4*pi * (1 - cos(R*Kn)) ./ k3;
%Green_2(1) = 0; % finite limit at k3==0
%Scaling_zero = -1i .* Green_2 .* Screening .* Window_m2;
%H0_2 = cell(3,1);
%H0_2{1} = Scaling_zero .* HH{2};
%H0_2{2} = -Scaling_zero .* HH{1};
%H0_2{3} = 0*HH{3};

%A = abs(Green_1);
%B = abs(Green_2);

%n = floor(numel(A)/2);

%clf
%subplot(2,1,1)
%loglog(k3(1:n), A(1:n), '^'), hold on
%loglog(k3(1:n), B(1:n), 'o')
%loglog(k3(1:n), 30./k3(1:n), 'k--')
%loglog(k3(1:n), 30./k3(1:n).^2, 'k:')
%xlim([1e0 2e2])
%ylim([1e-3 1e2])
%legend('|H_1|', '|H_2|', 'Location', 'best');
%title('H')

%A = abs(k3 .* Green_1);
%B = abs(Green_2);
%subplot(2,1,2)
%loglog(k3(1:n), A(1:n), '^'), hold on
%loglog(k3(1:n), B(1:n), 'o')
%loglog(k3(1:n), 30./k3(1:n), 'k--')
%xlim([1e0 2e2])
%ylim([1e-1 1e2])
%xlabel('Wavenumber (in free direction)');
%legend('|K_1H_1|', '|K_2H_2|', 'Location', 'best');
%title('KH')

% ------------------------------------------------------------------------------
function F = kaiser_exact_ft(k2, b2, w, scaling)

t = sqrt(b2 - k2*w^2);
F = 2*w*sinh(t)./t * scaling;
