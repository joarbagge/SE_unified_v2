function [G, Gl, G0] = SE1P_Stresslet_scaling(H, Hl, H0, opt, pre_fft)
%SE1P_Stresslet_scaling  Fourier-space scaling for 1-periodic stresslet
%
%   [G, Gl, G0] = SE1P_Stresslet_scaling(H, Hl, H0, opt, pre_fft)
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
k2 = SE_aft_k_vectors(opt.grid2, opt.box2, opt.actual_upsampling_global(1)); % free direction
k3 = SE_aft_k_vectors(opt.grid3, opt.box3, opt.actual_upsampling_global(2)); % free direction
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
    H0_cache{j} = squeeze(H{j}(1,:,:));
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
if numel(opt.local_modes1) > 0
  n1 = opt.local_modes1;
  k2 = SE_aft_k_vectors(opt.grid2, opt.box2, opt.actual_upsampling_local(1));
  k3 = SE_aft_k_vectors(opt.grid3, opt.box3, opt.actual_upsampling_local(2));
  [K{1}, K{2}, K{3}] = ndgrid(k1(n1), k2, k3);
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
      f2 = kaiser_exact_ft(k2.^2, b2, w, opt.kaiser_scaling);
      f3 = kaiser_exact_ft(k3.^2, b2, w, opt.kaiser_scaling);
      [F1, F2, F3] = ndgrid(f1(n1).^pw, f2.^pw, f3.^pw);
      F = F1.*F2.*F3; % tensor product of spatial directions
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

% Scaling for the zero mode [k1==0] (Vico method)
K{1} = 0;
k2 = SE_aft_k_vectors(opt.grid2, opt.box2, opt.actual_upsampling_zero(1));
k3 = SE_aft_k_vectors(opt.grid3, opt.box3, opt.actual_upsampling_zero(2));
[K{2}, K{3}] = ndgrid(k2, k3);
KK = K{2}.^2 + K{3}.^2;
C = KK/(4*opt.xi^2);
Kn = sqrt(KK);
R = opt.greens_truncation_R;
% Modified Green's function
%   - if opt.vico_variant == 1, base the stresslet completely on
%     the biharmonic
%   - if opt.vico_variant == 2, base the stresslet on a mix of
%     the harmonic and biharmonic (default)
% In any case, we will need the biharmonic
%lB = R; aB = 0; % half-optimal value for Vico formula convergence
%lB = R*exp(1/2); aB = -0.5*R^2; % optimal value for Vico formula convergence
%Biharmonic = 8*pi*( ...
%                    (besselj(0,R*Kn)-1)./(KK.*KK) ...
%                    + (R*(log(R/lB)+1)*besselj(1,R*Kn))./(KK.*Kn) ...
%                    - (R^2*(2*log(R/lB)+1)*besselj(0,R*Kn))./(4*KK) ...
%                    - (R*(R^2*log(R/lB)-aB)*besselj(1,R*Kn))./(4*Kn) ...
%                  );
%Biharmonic(1,1) = (pi*R^4/8) * (1-4*log(R/lB)) + pi*aB*R^2; % finite limit at k2==k3==0
% Use simplified code for the optimal case, for better speed
Biharmonic = 8*pi*( ...
                    (besselj(0,R*Kn)-1)./(KK.*KK) ...
                    + (0.5*R*besselj(1,R*Kn))./(KK.*Kn) ...
                  );
Biharmonic(1,1) = (-1/8)*pi*R^4; % finite limit at k2==k3==0
% We might also need the harmonic (if opt.vico_variant == 2)
if ~isfield(opt, 'vico_variant') || opt.vico_variant == 2
  %lH = R; % optimal value for Vico formula convergence
  %Harmonic = 4*pi*((1-besselj(0,R*Kn))./KK - R*log(R/lH)*besselj(1,R*Kn)./Kn);
  %Harmonic(1,1) = pi*R^2 * (1 - 2*log(R/lH)); % finite limit at k2==k3==0
  % Simplified code for better speed
  Harmonic = 4*pi*(1-besselj(0,R*Kn))./KK;
  Harmonic(1,1) = pi*R^2; % finite limit at k2==k3==0
end
if strcmp(opt.window, 'gaussian')
  Combo = (1+C) .* exp(-C*(1-eta*pw/2));
else
  Screening = (1+C) .* exp(-C);
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
  Combo = Screening .* Window_m2;
end
Scaling_zeroB = Biharmonic .* Combo;
if ~isfield(opt, 'vico_variant') || opt.vico_variant == 2
  Scaling_zeroH = Harmonic .* Combo;
end
if ~isempty(H0{1})
  H0_cache = H0;
end
% Relate harmonic/biharmonic to stresslet
G0_cache = cell(3,1);
HdotK = cell(3,1);
for j=1:3
  HdotK{j} = H0_cache{j,2}.*K{2} + H0_cache{j,3}.*K{3};
end
KdotHdotK = K{2}.*HdotK{2} + K{3}.*HdotK{3};
if ~isfield(opt, 'vico_variant') || opt.vico_variant == 2
  traceH = H0_cache{2,2} + H0_cache{3,3};
else
  traceH = H0_cache{1,1} + H0_cache{2,2} + H0_cache{3,3};
end
% Relation between stresslet and biharmonic
for j=1:3
  if j == 1 && (~isfield(opt, 'vico_variant') || opt.vico_variant == 2)
    continue;
  end
  KdotH = K{2}.*H0_cache{2,j} + K{3}.*H0_cache{3,j};
  vector_part = (HdotK{j} + KdotH + K{j}.*traceH).*KK - 2*K{j}.*KdotHdotK;
  G0_cache{j} = -1i * Scaling_zeroB .* vector_part;
end
% Relation between stresslet and harmonic
if ~isfield(opt, 'vico_variant') || opt.vico_variant == 2
  for j=1:3
    vector_part = K{j}.*H0_cache{1,1};
    if j == 1
      KdotH = K{2}.*H0_cache{2,j} + K{3}.*H0_cache{3,j};
      vector_part = vector_part + KdotH + HdotK{j};
      G0_cache{j} = 0;
    end
    G0_cache{j} = G0_cache{j} + 2i * Scaling_zeroH .* vector_part;
  end
end
if isempty(H0{1})
  % Store in G instead
  for j=1:numel(G)
    G{j}(1,:,:) = G0_cache{j};
  end
  G0 = {[], [], []}';
else
  G0 = G0_cache;
end

% ------------------------------------------------------------------------------
function F = kaiser_exact_ft(k2, b2, w, scaling)

t = sqrt(b2 - k2*w^2);
F = 2*w*sinh(t)./t * scaling;
