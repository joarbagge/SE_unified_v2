function [k1,k2,k3] = SE_k_vectors(M, box, shifted)
% SE_K_VECTORS  Compute k-vectors
%    [k1,k2,k3] = SE_K_VECTORS(M,BOX)
%    [k1,k2,k3] = SE_K_VECTORS(M,BOX,'shifted')

if nargin >= 3 && strcmp(shifted, 'shifted')
  shifted = true;
else
  shifted = false;
end

k1 = k_vec(M(1), box(1), shifted);
k2 = k_vec(M(2), box(2), shifted);
k3 = k_vec(M(3), box(3), shifted);

% ------------------------------------------------------------------------------
function k = k_vec(M, L, shifted)

if mod(M,2) == 0
  Mh = M/2;
  if ~shifted
    k = (2*pi/L) * [0:(Mh-1), -Mh:-1];
  else
    k = (2*pi/L) * [-Mh:(Mh-1)];
  end
else
  Mh = (M-1)/2;
  if ~shifted
    k = (2*pi/L) * [0:Mh, -Mh:-1];
  else
    k = (2*pi/L) * [-Mh:Mh];
  end
end
