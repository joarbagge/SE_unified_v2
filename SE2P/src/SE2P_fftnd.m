function [x, xl, x0] = SE2P_fftnd(x, opt)
%SE2P_fftnd  3-dimensional discrete Fourier Transform.
%
%   [x, xl, x0] = SE2P_fftnd(x, opt)
%      Returns the 3-dimensional DFT of the 3D array x with two periodic
%      directions and one free direction, in the arrays x, xl and
%      x0. The output xl contains the "local pad" result and x0
%      contains the zero mode.
%
%   Note: For simplicity, the local pad contains the zero mode as
%   well but it will be updated by the zero mode result later on.
%
%   Input parameters:
%       :param x: input vector of size (Mx,My,Mz)
%       :param opt.grid: vector containing the sizes [Mx,My,Mz]
%       :param opt.grid3: size of third dimension of x (i.e. Mz)
%       :param opt.actual_upsampling_zero: upsampling factor for the zero mode
%       :param opt.actual_upsampling_local: upsampling factor on the "local pad"
%       :param opt.actual_upsampling_global: upsampling factor for the whole domain
%       :param opt.local_modes1: list of modes to include in "local pad", x direction
%       :param opt.local_modes2: list of modes to include in "local pad", y direction
%
%   Return values:
%       :returns: **x** -- output vector of size (Mx,My,Mz*sg)
%       :returns: **xl** -- output vector of size (nx,ny,Mz*sl)
%       :returns: **x0** -- output vector of size (1,1,Mz*s0)
%
%   Here, nx=numel(opt.local_modes1), ny=numel(opt.local_modes2),
%   sg=opt.actual_upsampling_global, sl=opt.actual_upsampling_local,
%   s0=opt.actual_upsampling_zero.
%
%   The periodic directions are x and y, while the free direction is z.

sg = opt.actual_upsampling_global;
sl = opt.actual_upsampling_local;
s0 = opt.actual_upsampling_zero;

if sg == sl && sg == s0 % all upsampling factors are the same
  x = fftn(x, [opt.grid(1) opt.grid(2) round(opt.grid3*sg)]);
  xl = [];
  x0 = [];
else
  % Global domain
  Fxy = fft2(x, opt.grid(1), opt.grid(2));          % 2D fft in x and y
  x = fft(Fxy, round(opt.grid3*sg), 3);             % 1D fft in z (global)

  % Local pad
  L1 = opt.local_modes1;
  L2 = opt.local_modes2;
  xl = fft(Fxy(L1,L2,:), round(opt.grid3*sl), 3);   % 1D fft in z (local pad)

  % Zero mode
  F0 = Fxy(1,1,:); F0 = F0(:);
  x0 = fft(F0, round(opt.grid3*s0));                % 1D fft in z (zero mode)
end

end
