function pre_fft = SE1P_window_fft_precomp(opt)

% We assume periodicity in x and free in y and z.

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
Lx = opt.box(1);
x = 0:h:(Lx-h);
fx = h*window(x-Lx/2);
Fx = fft(fx).^pw;

% yz plane, whole domain, upsamping sg
Ly = opt.box2;
Lz = opt.box3;
My = opt.grid2;
Mz = opt.grid3;
sg = opt.actual_upsampling_global;
Myu = round(sg(1)*My);
Mzu = round(sg(2)*Mz);
y = 0:h:(Myu*h-h);
fy = h*window(y-Myu*h/2);
Fy = fft(fy, Myu).^pw;
z = 0:h:(Mzu*h-h);
fz = h*window(z-Mzu*h/2);
Fz = fft(fz, Mzu).^pw;
[F1, F2, F3] = ndgrid(Fx, Fy, Fz);
Fprod = F1.*F2.*F3; % tensor product of spatial directions
pre_fft.window = 1./Fprod;

% yz plane, local pad, upsampling sl
nx = opt.local_modes1;
sl = opt.actual_upsampling_local;
Myu = round(sl(1)*My);
Mzu = round(sl(2)*Mz);
y = 0:h:(Myu*h-h);
fy = h*window(y-Myu*h/2);
Fy = fft(fy, Myu).^pw;
z = 0:h:(Mzu*h-h);
fz = h*window(z-Mzu*h/2);
Fz = fft(fz, Mzu).^pw;
[F1, F2, F3] = ndgrid(Fx(nx), Fy, Fz);
Fprod = F1.*F2.*F3; % tensor product of spatial directions
pre_fft.window_local = 1./Fprod;

% yz plane, zero mode, upsampling s0
s0 = opt.actual_upsampling_zero;
Myu = round(s0(1)*My);
Mzu = round(s0(2)*Mz);
y = 0:h:(Myu*h-h);
fy = h*window(y-Myu*h/2);
Fy = fft(fy', Myu).^pw;
z = 0:h:(Mzu*h-h);
fz = h*window(z-Mzu*h/2);
Fz = fft(fz', Mzu).^pw;
F1 = Fx(1); [F2, F3] = ndgrid(Fy, Fz);
Fprod = F1.*F2.*F3; % tensor product of spatial directions
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
