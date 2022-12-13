function [G, Gl, G0] = SE2P_Stresslet_scaling(H, Hl, H0, opt, pre_fft)
%SE2P_Stresslet_scaling  Fourier-space scaling for 2-periodic stresslet
%
%   [G, Gl, G0] = SE2P_Stresslet_scaling(H, Hl, H0, opt, pre_fft)
%
%   Input parameters:
%       :param [H, Hl, H0]: input data in Fourier space (cell arrays)
%       :param opt: structure with Ewald options
%       :param pre_fft: optional precomputation structure:
%       :param pre_fft.window: W^(-pw) where W is the Fourier transform of the window function
%       :param pre_fft.window_local: W^(-pw) on the local pad
%       :param pre_fft.window_zero: W^(-pw) for the zero mode
%
%   :returns: **[G, Gl, G0]** -- output data in Fourier space (cell arrays)

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

% For the stresslet, we first use the biharmonic Green's function,
% and then use a relation between it and the stresslet.

% Scaling for the global Fourier domain
Biharmonic = -8*pi./(KK.*KK);
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
% Use relation between stresslet and biharmonic Green's function
G = cell(3,1);
HdotK = cell(3,1);
for j=1:3
  HdotK{j} = H{j,1}.*K{1} + H{j,2}.*K{2} + H{j,3}.*K{3};
end
KdotHdotK = K{1}.*HdotK{1} + K{2}.*HdotK{2} + K{3}.*HdotK{3};
for j=1:3
  KdotH = K{1}.*H{1,j} + K{2}.*H{2,j} + K{3}.*H{3,j};
  traceH = H{1,1} + H{2,2} + H{3,3};
  vector_part = -1i*(HdotK{j} + KdotH + K{j}.*traceH).*KK + 2i*K{j}.*KdotHdotK;
  G{j} = Scaling .* vector_part;
end

% Scaling for the local pad (same structure as above)
if numel(opt.local_modes1) > 0 || numel(opt.local_modes2) > 0
  n1 = opt.local_modes1;
  n2 = opt.local_modes2;
  k3 = SE_aft_k_vectors(opt.grid3, opt.box3, opt.actual_upsampling_local);
  [K{1}, K{2}, K{3}] = ndgrid(k1(n1), k2(n2), k3);
  KK = K{1}.^2 + K{2}.^2 + K{3}.^2;
  Biharmonic = -8*pi./(KK.*KK);
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
  % Use relation between stresslet and biharmonic Green's function
  Gl = cell(3,1);
  HdotK = cell(3,1);
  for j=1:3
    HdotK{j} = Hl{j,1}.*K{1} + Hl{j,2}.*K{2} + Hl{j,3}.*K{3};
  end
  KdotHdotK = K{1}.*HdotK{1} + K{2}.*HdotK{2} + K{3}.*HdotK{3};
  for j=1:3
    KdotH = K{1}.*Hl{1,j} + K{2}.*Hl{2,j} + K{3}.*Hl{3,j};
    traceH = Hl{1,1} + Hl{2,2} + Hl{3,3};
    vector_part = -1i*(HdotK{j} + KdotH + K{j}.*traceH).*KK + 2i*K{j}.*KdotHdotK;
    Gl{j} = Scaling_local .* vector_part;
  end
else
  Gl = {[], [], []}';
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
% Modified Green's function (based on -8*pi/k3)
if ~isfield(opt, 'vico_variant') || opt.vico_variant == 4
  Green = -8*pi * (1 - cos(R*Kn)) ./ k3;
  Green(1) = 0; % finite limit at k3==0
  Kfactor = -1i;
elseif opt.vico_variant == 3
  Green = -8*pi * (1 - cos(R*Kn) - R*Kn.*sin(R*Kn)) ./ KK;
  Green(1) = 4*pi*R^2; % finite limit at k3==0
  Kfactor = -1i*k3;
elseif opt.vico_variant == 2
  Green = -4*pi*( ...
                  -2*R*Kn.*sin(R*Kn) ...
                  + (R^2*KK-2).*cos(R*Kn) ...
                  + 2 ...
                ) ./ (KK.*k3);
  Green(1) = 0; % finite limit at k3==0
  Kfactor = -1i*KK;
elseif opt.vico_variant == 1
  Green = -4*pi*( ...
                  R*Kn.*(R^2*KK-6).*sin(R*Kn) ...
                  + 3*(R^2*KK-2).*cos(R*Kn) ...
                  + 6 ...
                ) ./ (3*KK.*KK);
  Green(1) = -pi*R^4/3; % finite limit at k3==0
  Kfactor = -1i*k3.*KK;
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
% Use simplified relation between stresslet and biharmonic Green's function
G0_cache = cell(3,1);
G0_cache{1} = Scaling_zero .* (H0_cache{1,3} + H0_cache{3,1});
G0_cache{2} = Scaling_zero .* (H0_cache{2,3} + H0_cache{3,2});
G0_cache{3} = Scaling_zero .* (H0_cache{1,1} + H0_cache{2,2} + H0_cache{3,3});
if isempty(H0{1})
  % Store in G instead
  for j=1:numel(G)
    G{j}(1,1,:) = G0_cache{j};
  end
  G0 = {[], [], []}';
else
  G0 = G0_cache;
end

% (1), regularizing -8*pi/k^4 (whole biharmonic)
%Green_1 = -4*pi*( ...
%               R*Kn.*(R^2*KK-6).*sin(R*Kn) ...
%               + 3*(R^2*KK-2).*cos(R*Kn) ...
%               + 6 ...
%             ) ./ (3*KK.*KK);
%Green_1(1) = -pi*R^4/3; % finite limit at k3==0
%Kfac_1 = -1i*k3.*KK; % remaining factor of K_jlm
%Scaling_zero = Kfac_1 .* Green_1 .* Screening .* Window_m2;
%G0_1 = cell(3,1);
%G0_1{1} = Scaling_zero .* (H0_cache{1,3} + H0_cache{3,1});
%G0_1{2} = Scaling_zero .* (H0_cache{2,3} + H0_cache{3,2});
%G0_1{3} = Scaling_zero .* (H0_cache{1,1} + H0_cache{2,2} + H0_cache{3,3});
% G0_1 and G0 agree!

% (2) regularizing -8*pi/k^3
%Green_2 = -4*pi*( ...
%               -2*R*Kn.*sin(R*Kn) ...
%               + (R^2*KK-2).*cos(R*Kn) ...
%               + 2 ...
%             ) ./ (KK.*k3);
%Green_2(1) = 0; % finite limit at k3==0
%Kfac_2 = -1i*KK; % remaining factor of K_jlm
%Scaling_zero = Kfac_2 .* Green_2 .* Screening .* Window_m2;
%G0_2 = cell(3,1);
%G0_2{1} = Scaling_zero .* (H0_cache{1,3} + H0_cache{3,1});
%G0_2{2} = Scaling_zero .* (H0_cache{2,3} + H0_cache{3,2});
%G0_2{3} = Scaling_zero .* (H0_cache{1,1} + H0_cache{2,2} + H0_cache{3,3});
% G0_2 also agrees relatively well with the other ones!

% (3) regularizing -8*pi/k^2
%Green_3 = -8*pi * (1 - cos(R*Kn) - R*Kn.*sin(R*Kn)) ./ KK;
%Green_3(1) = 4*pi*R^2; % finite limit at k3==0
%Kfac_3 = -1i*k3; % remaining factor of K_jlm
%Scaling_zero = Kfac_3 .* Green_3 .* Screening .* Window_m2;
%G0_3 = cell(3,1);
%G0_3{1} = Scaling_zero .* (H0_cache{1,3} + H0_cache{3,1});
%G0_3{2} = Scaling_zero .* (H0_cache{2,3} + H0_cache{3,2});
%G0_3{3} = Scaling_zero .* (H0_cache{1,1} + H0_cache{2,2} + H0_cache{3,3});
% G0_3 does NOT agree with G0 and the previous G0_1 and G0_2!

% (4) regularizing -8*pi/k
%Green_4 = -8*pi * (1 - cos(R*Kn)) ./ k3;
%Green_4(1) = 0; % finite limit at k3==0
%Kfac_4 = -1i; % remaining factor of K_jlm
%Scaling_zero = Kfac_4 .* Green_4 .* Screening .* Window_m2;
%G0_4 = cell(3,1);
%G0_4{1} = Scaling_zero .* (H0_cache{1,3} + H0_cache{3,1});
%G0_4{2} = Scaling_zero .* (H0_cache{2,3} + H0_cache{3,2});
%G0_4{3} = Scaling_zero .* (H0_cache{1,1} + H0_cache{2,2} + H0_cache{3,3});
% G0_4 does NOT agree with G0, but it *does* agree with G0_3

%A = abs(Green_1);
%B = abs(Green_2);
%C = abs(Green_3);
%D = abs(Green_4);

%n = floor(numel(A)/2);

%clf
%subplot(2,1,1)
%loglog(k3(1:n), A(1:n), '.-'), hold on
%loglog(k3(1:n), B(1:n), '.-')
%loglog(k3(1:n), C(1:n), '^')
%loglog(k3(1:n), D(1:n), 'o')
%loglog(k3(1:n), 30./k3(1:n), 'k--')
%loglog(k3(1:n), 30./k3(1:n).^2, 'k:')
%xlim([1e0 2e2])
%ylim([1e-3 1e2])
%legend('|B_1|', '|B_2|', '|B_3|', '|B_4|', 'Location', 'best');
%title('B')

%A = abs(Kfac_1 .* Green_1);
%B = abs(Kfac_2 .* Green_2);
%C = abs(Kfac_3 .* Green_3);
%D = abs(Kfac_4 .* Green_4);
%subplot(2,1,2)
%loglog(k3(1:n), A(1:n), '.-'), hold on
%loglog(k3(1:n), B(1:n), '.-')
%loglog(k3(1:n), C(1:n), '^')
%loglog(k3(1:n), D(1:n), 'o')
%loglog(k3(1:n), 30*k3(1:n), 'k--')
%loglog(k3(1:n), 30./k3(1:n), 'k:')
%xlim([1e0 2e2])
%ylim([1e-1 1e4])
%xlabel('Wavenumber (in free direction)');
%legend('|K_1B_1|', '|K_2B_2|', '|K_3B_3|', '|K_4B_4|', 'Location', 'best');
%title('KB')

% ------------------------------------------------------------------------------
function F = kaiser_exact_ft(k2, b2, w, scaling)

t = sqrt(b2 - k2*w^2);
F = 2*w*sinh(t)./t * scaling;
