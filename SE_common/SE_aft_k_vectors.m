function k = SE_aft_k_vectors(M, L, upsampling)
%SE_aft_k_vectors  Compute k-vectors for the adaptive Fourier grid
%
%   k = SE_aft_k_vectors(M, box, upsampling)
%
%   Input parameters:
%       :param M: grid size in a particular direction (scalar)
%       :param L: box side length in same direction (scalar)
%       :param upsampling: upsampling factor (1 = no upsampling)
%
%   :returns: **k** -- k-vectors in given direction

  Mu = round(M*upsampling);

  if mod(Mu,2) == 0
      kint = (-Mu/2):(Mu/2-1);
  else
      kint = -(Mu-1)/2:(Mu-1)/2;
  end

  kint = fftshift(kint); % standard reordering, so that kint(1)==0

  h = L/M;
  Lu = Mu*h;
  k = 2*pi*kint/Lu;
  if (abs(k(1)) > eps)
      k = circshift(k, [1 1]);
  end
  assert(abs(k(1)) < eps);

end
