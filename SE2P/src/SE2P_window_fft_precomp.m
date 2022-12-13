function pre_fft = SE2P_window_fft_precomp(opt)

% We assume periodicity in x and y and free in z.

pw = opt.window_scaling_power; % should be 2 or 1 typically
sh = opt.window_shape;
w = opt.window_halfwidth;
h = opt.h;

if strcmp(opt.window, 'gaussian')
  error('FFT precomputing not supported for Gaussian window');
elseif strcmp(opt.window, 'expsemicirc')
  window = @(x) expsemicirc(x, sh, w);
elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
  window = @(x) kaiser_exact(x, sh, w, opt.kaiser_scaling);
else
  error('Unsupported window function');
end

% Common in all parts
F = cell(1,2);
for d=1:2
  L = opt.box(d);
  x = 0:h:(L-h);
  f = h*window(x-L/2);
  F{d} = fft(f).^pw;
end

% z direction, whole domain, upsampling sg
Lz = opt.box3;
Mz = opt.grid3;
sg = opt.actual_upsampling_global;
Mzu = round(sg*Mz);
z = 0:h:(Mzu*h-h);
fz = h*window(z-Mzu*h/2);
Fz = fft(fz, Mzu).^pw;
[F1, F2, F3] = ndgrid(F{1}, F{2}, Fz);
Fprod = F1.*F2.*F3; % tensor product of spatial directions
pre_fft.window = 1./Fprod;

% z direction, local pad, upsampling sl
nx = opt.local_modes1;
ny = opt.local_modes2;
sl = opt.actual_upsampling_local;
Mzu = round(sl*Mz);
z = 0:h:(Mzu*h-h);
fz = h*window(z-Mzu*h/2);
Fz = fft(fz, Mzu).^pw;
[F1, F2, F3] = ndgrid(F{1}(nx), F{2}(ny), Fz);
Fprod = F1.*F2.*F3; % tensor product of spatial directions
pre_fft.window_local = 1./Fprod;

% z direction, zero mode, upsampling s0
s0 = opt.actual_upsampling_zero;
Mzu = round(s0*Mz);
z = 0:h:(Mzu*h-h);
fz = h*window(z-Mzu*h/2);
Fz = fft(fz', Mzu).^pw;
F1 = F{1}(1); F2 = F{2}(1); F3 = Fz;
Fprod = F1*F2*F3; % tensor product of spatial directions
pre_fft.window_zero = 1./Fprod;

end

% ------------------------------------------------------------------------------
function z = expsemicirc(x, beta, w)

t = sqrt(1-(x/w).^2);
z = exp(beta*(t-1)) .* (abs(x) <= w);

end

% ------------------------------------------------------------------------------
function z = kaiser_exact(x, beta, w, scaling)

t = sqrt(1-(x/w).^2);
I = besseli(0,beta*t) * scaling;
z = I .* (abs(x) <= w);

end
