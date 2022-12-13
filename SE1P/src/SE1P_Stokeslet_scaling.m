function [H, Hl, H0] = SE1P_Stokeslet_scaling(H, Hl, H0, opt, pre_fft)
%SE1P_Stokeslet_scaling  Fourier-space scaling for 1-periodic stokeslet
%
%   [H, Hl, H0] = SE1P_Stokeslet_scaling(H, Hl, H0, opt, pre_fft)
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
    H0_cache{j} = squeeze(H{j}(1,:,:));
  end
end
% Use relation between stokeslet and biharmonic Green's function
KdotH = K1.*H{1} + K2.*H{2} + K3.*H{3};
H{1} = Scaling .* (KK.*H{1} - KdotH.*K1);
H{2} = Scaling .* (KK.*H{2} - KdotH.*K2);
H{3} = Scaling .* (KK.*H{3} - KdotH.*K3);

% Scaling for the local pad (same structure as above)
if numel(opt.local_modes1) > 0
  n1 = opt.local_modes1;
  k2 = SE_aft_k_vectors(opt.grid2, opt.box2, opt.actual_upsampling_local(1));
  k3 = SE_aft_k_vectors(opt.grid3, opt.box3, opt.actual_upsampling_local(2));
  [K1, K2, K3] = ndgrid(k1(n1), k2, k3);
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
      f2 = kaiser_exact_ft(k2.^2, b2, w, opt.kaiser_scaling);
      f3 = kaiser_exact_ft(k3.^2, b2, w, opt.kaiser_scaling);
      [F1, F2, F3] = ndgrid(f1(n1).^pw, f2.^pw, f3.^pw);
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

% Scaling for the zero mode [k1==0] (Vico method)
k2 = SE_aft_k_vectors(opt.grid2, opt.box2, opt.actual_upsampling_zero(1));
k3 = SE_aft_k_vectors(opt.grid3, opt.box3, opt.actual_upsampling_zero(2));
[K2, K3] = ndgrid(k2, k3);
KK = K2.^2 + K3.^2;
C = KK/(4*opt.xi^2);
Kn = sqrt(KK);
R = opt.greens_truncation_R;
% Modified Green's function
%   - if opt.vico_variant == 1, base the stokeslet completely on
%     the biharmonic
%   - if opt.vico_variant == 2, base the stokeslet on a mix of
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
% Relate harmonic/biharmonic to stokeslet
SKK = -KK .* Scaling_zeroB;
% Use harmonic or biharmonic for the first component
if ~isfield(opt, 'vico_variant') || opt.vico_variant == 2
  H0_cache{1} = 2 * Scaling_zeroH .* H0_cache{1};
elseif opt.vico_variant == 1
  H0_cache{1} = SKK .* H0_cache{1};
end
% Always use biharmonic for the other two components
KdotH = K2.*H0_cache{2} + K3.*H0_cache{3};
Scaling_zeroB = Scaling_zeroB .* KdotH;
H0_cache{2} = SKK .* H0_cache{2} + Scaling_zeroB .* K2;
H0_cache{3} = SKK .* H0_cache{3} + Scaling_zeroB .* K3;
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
