function [x, xl, x0] = SE1P_fftnd(x, opt)
%SE1P_fftnd  3-dimensional discrete Fourier Transform.
%
%   [x, xl, x0] = SE1P_fftnd(x, opt)
%      Returns the 3-dimensional DFT of the 3D array x with one periodic
%      directions and two free direction, in the arrays x, xl and
%      x0. The output xl contains the "local pad" result and x0
%      contains the zero mode.
%
%   Input parameters:
%       :param x: input vector of size (Mx,My,Mz)
%       :param opt.grid: vector containing the sizes [Mx,My,Mz]
%       :param opt.grid2: size of second dimension of x (i.e. My)
%       :param opt.grid3: size of third dimension of x (i.e. Mz)
%       :param opt.actual_upsampling_zero: upsampling factor for the zero mode
%       :param opt.actual_upsampling_local: upsampling factor on the "local pad"
%       :param opt.actual_upsampling_global: upsampling factor for the whole domain
%       :param opt.local_modes1: list of modes to include in "local pad"
%
%   Return values:
%       :returns: **x** -- output vector of size (Mx,My*sg(1),Mz*sg(2))
%       :returns: **xl** -- output vector of size (nx,My*sl(1),Mz*sl(2))
%       :returns: **x0** -- output vector of size (1,My*s0(1),Mz*s0(2))
%
%   Here, nx=numel(opt.local_modes1), sg=opt.actual_upsampling_global,
%   sl=opt.actual_upsampling_local, s0=opt.actual_upsampling_zero.
%
%   The periodic direction is x, while the free directions are y and z.

sg = opt.actual_upsampling_global;
sl = opt.actual_upsampling_local;
s0 = opt.actual_upsampling_zero;

if all(sg == sl) && all(sg == s0) % all upsampling factors are the same
  x = fftn(x, [opt.grid(1) round(opt.grid2*sg(1)) round(opt.grid3*sg(2))]);
  xl = [];
  x0 = [];
else
  % Global domain
  Fx = fft(x, opt.grid(1));                                % 1D fft in x
  x = fft(fft(Fx, round(opt.grid2*sg(1)), 2), ...
          round(opt.grid3*sg(2)), 3);                      % 2D fft in y and z (global)

  % Local pad
  L1 = opt.local_modes1;
  xl = fft(fft(Fx(L1,:,:), round(opt.grid2*sl(1)), 2), ...
           round(opt.grid3*sl(2)), 3);                     % 2D fft in y and z (local pad)

  % Zero mode
  F0 = squeeze(Fx(1,:,:));
  x0 = fft2(F0, round(opt.grid2*s0(1)), round(opt.grid3*s0(2))); % 2D fft in y and z (zero mode)
end

end
