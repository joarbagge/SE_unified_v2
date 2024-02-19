function pot = SE_real_space(kernel, varargin)
% Compute real-space part of the Ewald potential. The kernel may
% optionally be regularized (smoothened) according to Tlupova & Beale (2019),
% or Tlupova & Beale (2022).
%
% pot = SE_real_space('stokes_sl', opt, targets, sources, density)
% pot = SE_real_space('stokes_dl', opt, targets, sources, normals, density)
% pot = SE_real_space(..., delta)
% pot = SE_real_space('stokes_dl', ..., smoothing)
%
% The first argument is 'stokes_sl' for the Stokes single-layer
% potential (stokeslet), or 'stokes_dl' for the Stokes
% double-layer potential (stresslet).
%
% Input parameters:
% - opt: option struct, see further below
% - targets: Nt×3 matrix with target point locations, i.e., points
%   to evaluate the potential at
% - sources: Ns×3 matrix with source point locations
% - normals: (only for 'stokes_dl') Ns×3 matrix with surface unit
%   normals at source points
% - density: Ns×3 matrix with layer density values at source
%   points
% - delta: Ns×1 vector or scalar with regularization (smoothening)
%   distances (for a scalar, the same value is used for all sources).
%   If zero (the default), no regularization is performed, except
%   that the "self-interaction" is ignored.
% - smoothing: (optional, only for 'stokes_dl') which smoothing
%   function to use; one of 'auto' (the default), 'full' and
%   'simple'
%
% Valid options for the opt struct are as follows. All options
% are required.
% - opt.periodicity: the number of periodic directions (0, 1, 2 or 3)
% - opt.box: size of primary cell [L1,L2,L3], containing all
%   sources and targets
% - opt.xi: Ewald decomposition parameter
% - opt.rc: cutoff radius
%
% Periodicity is assumed in the first `opt.periodicity` directions,
% while the remaining `3-opt.periodicity` directions are assumed
% free.
%
% See also SE_Stokeslet_params and SE_Stresslet_params.
%
% Output:
% - pot: Nt×3 matrix with potential values at targets

% Note: This function is the Matlab reference implementation.
% It is relatively slow and O(N^2) since the distance between
% every target and every source is computed.
% See also: ../mex/SE_real_space_mex.c

%%
% Argument handling
double_layer = false;
if strcmp(kernel, 'stokes_sl')
    kernel_fun = @stokes_sl_kernel;
elseif strcmp(kernel, 'stokes_dl')
    double_layer = true;
    kernel_fun = @stokes_dl_kernel;
else
    error('Unsupported kernel %s', kernel);
end
if nargin < 5 + double_layer
    error('SE_real_space expects at least %d arguments, got %d', ...
          5+double_layer, nargin);
end
opt = varargin{1};
targets = varargin{2};
sources = varargin{3};
if double_layer
    normals = varargin{4};
else
    normals = [];
end
density = varargin{4+double_layer};
% Default values
delta = 0;
smoothing = 'auto';
% Optional arguments
if nargin >= 6 + double_layer
    delta = varargin{5+double_layer};
end
if nargin >= 7 + double_layer
    smoothing = varargin{6+double_layer};
end
if ~double_layer && ~strcmp(smoothing, 'auto')
    error('%s does not support smoothing "%s"', kernel, smoothing);
end
% End argument handling

%%
% Select smoothing function
num_terms = 4;
if ~double_layer
    smooth(1).fun = @smoothfun1;
    % For smoothfun1, 1 Taylor term is enough to get machine precision
    [smooth(1).taylor, smooth(1).threshold] = smoothfun1_taylor_setup(1);
    smooth(2).fun = @smoothfun2;
    [smooth(2).taylor, smooth(2).threshold] = smoothfun2_taylor_setup(num_terms);
else
    if strcmp(smoothing, 'auto')
        % Density subtraction is not used; must use the full
        % smoothing function from Tlupova & Beale (2019).
        smoothing = 'full';
    end
    if strcmp(smoothing, 'simple')
        smooth.fun = @smoothfun3_simple;
        [smooth.taylor, smooth.threshold] = smoothfun3_simple_taylor_setup(num_terms);
    elseif strcmp(smoothing, 'full')
        smooth.fun = @smoothfun3_full;
        [smooth.taylor, smooth.threshold] = smoothfun3_full_taylor_setup(num_terms);
    elseif double_layer
        error('Unknown smoothing function "%s"', smoothing);
    end
end

%%
% Check opt struct
assert(isfield(opt, 'periodicity'), ['Number of periodic directions "periodicity"' ...
                                     ' must be given in opt struct']);
assert(isfield(opt, 'box'), 'Primary cell size "box" must be given in opt struct');
assert(isfield(opt, 'xi'), 'Ewald decomposition parameter "xi" must be given in opt struct');
assert(isfield(opt, 'rc'), 'Cutoff radius "rc" must be given in opt struct');

%%
% Call functions for each periodic case
Nt = size(targets, 1);
pot = zeros(Nt, 3);
if opt.periodicity == 0
    pot = core_0p(pot, kernel_fun, opt, targets, sources, normals, density, delta, smooth);
elseif opt.periodicity == 1
    pot = core_1p(pot, kernel_fun, opt, targets, sources, normals, density, delta, smooth);
elseif opt.periodicity == 2
    pot = core_2p(pot, kernel_fun, opt, targets, sources, normals, density, delta, smooth);
elseif opt.periodicity == 3
    pot = core_3p(pot, kernel_fun, opt, targets, sources, normals, density, delta, smooth);
else
    error('Number of periodic directions must be in {0,1,2,3}, found %d', opt.periodicity);
end

end % function SE_real_space

function pot = core_0p(pot, kernel_fun, opt, targets, sources, normals, density, delta, smooth)
    pot = core_single_image(pot, kernel_fun, opt, targets, sources, ...
                            normals, density, delta, smooth);
end

function pot = core_1p(pot, kernel_fun, opt, targets, sources, normals, density, delta, smooth)
    for shift1=[-1,0,1]
        shift = [shift1,0,0] .* opt.box;
        trg = targets + shift;
        pot = core_single_image(pot, kernel_fun, opt, trg, sources, ...
                                normals, density, delta, smooth);
    end
end

function pot = core_2p(pot, kernel_fun, opt, targets, sources, normals, density, delta, smooth)
    for shift1=[-1,0,1]
        for shift2=[-1,0,1]
            shift = [shift1,shift2,0] .* opt.box;
            trg = targets + shift;
            pot = core_single_image(pot, kernel_fun, opt, trg, sources, ...
                                    normals, density, delta, smooth);
        end
    end
end

function pot = core_3p(pot, kernel_fun, opt, targets, sources, normals, density, delta, smooth)
    for shift1=[-1,0,1]
        for shift2=[-1,0,1]
            for shift3=[-1,0,1]
                shift = [shift1,shift2,shift3] .* opt.box;
                trg = targets + shift;
                pot = core_single_image(pot, kernel_fun, opt, trg, sources, ...
                                        normals, density, delta, smooth);
            end
        end
    end
end

function pot = core_single_image(pot, kernel_fun, opt, targets, sources, ...
                                 normals, density, delta, smooth)
    xi = opt.xi;
    xi2 = xi*xi;
    rc = opt.rc;
    rc2 = rc*rc;
    Nt = size(targets, 1);
    Ns = size(sources, 1);
    % Loop over targets
    for j=1:Nt
        target = targets(j,:);
        % Find sources within rc of target
        rvec = bsxfun(@minus, target, sources);
        r2 = sum(rvec.^2, 2);
        within_rc = (r2 <= rc2);
        if ~any(within_rc)
            continue;
        end
        % Slice data
        rvec = rvec(within_rc,:);
        r2 = r2(within_rc,:);
        if isempty(normals)
            n = [];
        else
            n = normals(within_rc,:);
        end
        f = density(within_rc,:);
        if numel(delta) == 1
            d = delta;
        else
            d = delta(within_rc);
        end
        % Common computations (regardless of kernel)
        r = sqrt(r2);
        xir = xi*r;
        xir2 = xi2*r2;
        xiexp = xi * exp(-xir2);
        rdotf = sum(rvec.*f, 2);
        A = erf(xir)./r;
        % Replace A by its limit close to zero (relative error at 1e-14
        % should be around 4e-29, so this is very accurate)
        S = sqrt(pi);
        A(xir < 1e-14) = 2*xi/S;
        B = 2*xiexp/S;
        C = (A-B)./r2;
        % Replace C by its Taylor expansion close to zero
        % (relative error at most around 3e-13 for four terms)
        % NOTE: This might be a bit unnecessary, as C goes into terms
        % that will anyway go to zero as r -> 0. But we do it
        % anyway for good measure.
        term1 = xi2*xi/S;
        %C(xir < 1.83e-4) = (4/3)*term1; % one term
        mask = (xir < 4.75e-2);
        term2 = xir2(mask).*term1;
        term3 = xir2(mask).*term2;
        term4 = xir2(mask).*term3;
        C(mask) = (4/3)*term1 - (4/5)*term2 + (2/7)*term3 - (2/27)*term4; % four terms
        % Kernel-dependent computations
        pot(j,:) = pot(j,:) + kernel_fun(xi, xi2, rvec, r, r2, n, f, rdotf, d, A, B, C, smooth);
    end
end

function pot = stokes_sl_kernel(xi, xi2, rvec, r, r2, n, f, rdotf, d, A, B, C, smooth)
    % Singular part, i.e., terms coming from the full (free-space) stokeslet
    if numel(d) == 1 && d == 0
        % Punctured trapezoidal rule
        s1 = 1 ./ r;
        s1(r < 1e-14) = 0; % remove point where r==0 (1e-14 is ad hoc)
        s2 = 1 ./ (r2.*r);
        s2(r < 1e-14) = 0; % remove point where r==0 (1e-14 is ad hoc)
    else
        % Beale regularization
        rd = r ./ d;
        s1 = smooth(1).fun(rd) ./ r;
        % Replace by known Taylor expansion for r close to zero
        I1 = (rd < smooth(1).threshold);
        if numel(d) == 1
            s1(I1) = smooth(1).taylor(d, rd(I1));
        else
            s1(I1) = smooth(1).taylor(d(I1), rd(I1));
        end
        s2 = smooth(2).fun(rd) ./ (r2.*r);
        % Replace by known Taylor expansion for r close to zero
        I2 = (rd < smooth(2).threshold);
        if numel(d) == 1
            s2(I2) = smooth(2).taylor(d, rd(I2));
        else
            s2(I2) = smooth(2).taylor(d(I2), rd(I2));
        end
    end
    t1 = bsxfun(@times, s1, f);
    t2 = bsxfun(@times, rdotf.*s2, rvec);
    % Non-singular part, i.e., terms coming from the Fourier-space part of the stokeslet
    tmp = -C.*rdotf;
    t3 = bsxfun(@times, tmp, rvec);
    t4 = bsxfun(@times, -B-A, f);
    % Sum up all terms
    pot = sum(t1+t2+t3+t4, 1) / (8*pi);
end

function pot = stokes_dl_kernel(xi, xi2, rvec, r, r2, n, f, rdotf, d, A, B, C, smooth)
    r4 = r2.*r2;
    rdotn = sum(rvec.*n, 2);
    fdotn = sum(f.*n, 2);
    % Singular part, i.e., terms coming from the full (free-space) stresslet
    if numel(d) == 1 && d == 0
        % Punctured trapezoidal rule
        s = 1 ./ (r4.*r);
        s(r < 1e-14) = 0; % remove point where r==0 (1e-14 is ad hoc)
    else
        % Beale regularization
        rd = r ./ d;
        s = smooth.fun(rd) ./ (r4.*r);
        % Replace by known Taylor expansion for r close to zero
        I = (rd < smooth.threshold);
        if numel(d) == 1
            s(I) = smooth.taylor(d, rd(I));
        else
            s(I) = smooth.taylor(d(I), rd(I));
        end
    end
    t1 = bsxfun(@times, -6*s.*rdotf.*rdotn, rvec);
    % Non-singular part, i.e., terms coming from the Fourier-space part of the stresslet
    D = 2*xi2*B;
    tmp = (6*C-2*D).*rdotf.*rdotn./r2;
    % At r==0, the term will become zero (1e-14 is ad hoc)
    tmp(r < 1e-14) = 0;
    t2 = bsxfun(@times, tmp + D.*fdotn, rvec);
    t3 = bsxfun(@times, D.*rdotn, f);
    t4 = bsxfun(@times, D.*rdotf, n);
    % Sum up all terms
    pot = sum(t1+t2+t3+t4, 1) / (8*pi);
end

function s1 = smoothfun1(x)
% Function s1# from Tlupova & Beale (2019) for use on-surface.
    a = erf(x);
    x2 = x.^2;
    b = 2*x2 - 5;
    c = -(2/3/sqrt(pi)) * x .* b .* exp(-x2);
    s1 = a + c;
end

function s1 = smoothfun1_taylor(d, rd, terms)
    term1 = 1./(sqrt(pi)*d);
    s1 = (16/3)*term1;
    if terms == 1
        return;
    end
    rd2 = rd.*rd;
    term2 = term1 .* rd2;
    s1 = s1 - (16/3)*term2;
    if terms == 2
        return;
    end
    term3 = term2 .* rd2;
    s1 = s1 + (16/5)*term3;
    if terms == 3
        return;
    end
    term4 = term3 .* rd2;
    s1 = s1 - (80/63)*term4;
    if terms == 4
        return;
    end
    error('%d terms not supported\n', terms);
end

function [fun, thres] = smoothfun1_taylor_setup(terms)
    fun = @(d, rd) smoothfun1_taylor(d, rd, terms);
    if terms == 1
        thres = 2.11e-8;
    elseif terms == 2
        thres = 1.65e-4;
    elseif terms == 3
        thres = 3.51e-3;
    elseif terms == 4
        thres = 1.68e-2;
    else
        error('%d terms not supported\n', terms);
    end
end

function s2 = smoothfun2(x)
% Function s2# from Tlupova & Beale (2019) for use on-surface.
    a = erf(x);
    x2 = x.^2;
    % Remember: x.^4 is SLOW in Matlab.
    b = 4*x2.*x2 - 14*x2 + 3;
    c = -(2/3/sqrt(pi)) * x .* b .* exp(-x2);
    s2 = a + c;
end

function s2 = smoothfun2_taylor(d, rd, terms)
    term1 = 1./(sqrt(pi)*d.^3);
    s2 = (32/3)*term1;
    if terms == 1
        return;
    end
    rd2 = rd.*rd;
    term2 = term1 .* rd2;
    s2 = s2 - (64/5)*term2;
    if terms == 2
        return;
    end
    term3 = term2 .* rd2;
    s2 = s2 + (160/21)*term3;
    if terms == 3
        return;
    end
    term4 = term3 .* rd2;
    s2 = s2 - (80/27)*term4;
    if terms == 4
        return;
    end
    error('%d terms not supported\n', terms);
end

function [fun, thres] = smoothfun2_taylor_setup(terms)
    fun = @(d, rd) smoothfun2_taylor(d, rd, terms);
    if terms == 1
        thres = 9.13e-5;
    elseif terms == 2
        thres = 2.21e-3;
    elseif terms == 3
        thres = 1.15e-2;
    elseif terms == 4
        thres = 3.18e-2;
    else
        error('%d terms not supported\n', terms);
    end
end

function s3 = smoothfun3_full(x)
% Function s3# from Tlupova & Beale (2019) for use on-surface.
    a = erf(x);
    x2 = x.^2;
    % Remember: x.^4 is SLOW in Matlab.
    x4 = x2.*x2;
    b = 8*x4.*x2 - 36*x4 + 6*x2 + 9;
    c = -(2/9/sqrt(pi)) * x .* b .* exp(-x2);
    s3 = a + c;
end

function s3 = smoothfun3_full_taylor(d, rd, terms)
    term1 = 1./(sqrt(pi)*d.^5);
    s3 = (128/15)*term1;
    if terms == 1
        return;
    end
    rd2 = rd.*rd;
    term2 = term1 .* rd2;
    s3 = s3 - (640/63)*term2;
    if terms == 2
        return;
    end
    term3 = term2 .* rd2;
    s3 = s3 + (160/27)*term3;
    if terms == 3
        return;
    end
    term4 = term3 .* rd2;
    s3 = s3 - (224/99)*term4;
    if terms == 4
        return;
    end
    error('%d terms not supported\n', terms);
end

function [fun, thres] = smoothfun3_full_taylor_setup(terms)
    fun = @(d, rd) smoothfun3_full_taylor(d, rd, terms);
    if terms == 1
        thres = 2.11e-3;
    elseif terms == 2
        thres = 1.05e-2;
    elseif terms == 3
        thres = 2.88e-2;
    elseif terms == 4
        thres = 5.78e-2;
    else
        error('%d terms not supported\n', terms);
    end
end

function s3 = smoothfun3_simple(x)
% Function s3# from Tlupova & Beale (2022) for use on-surface,
% together with density subtraction.
    a = erf(x);
    x2 = x.^2;
    % Remember: x.^4 is SLOW in Matlab.
    x4 = x2.*x2;
    b = -4*x4 + 6*x2 + 9;
    c = -(2/9/sqrt(pi)) * x .* b .* exp(-x2);
    s3 = a + c;
end

function s3 = smoothfun3_simple_taylor(d, rd, terms)
    term1 = 1./(sqrt(pi)*d.^5);
    s3 = (64/45)*term1;
    if terms == 1
        return;
    end
    rd2 = rd.*rd;
    term2 = term1 .* rd2;
    s3 = s3 - (80/63)*term2;
    if terms == 2
        return;
    end
    term3 = term2 .* rd2;
    s3 = s3 + (16/27)*term3;
    if terms == 3
        return;
    end
    term4 = term3 .* rd2;
    s3 = s3 - (56/297)*term4;
    if terms == 4
        return;
    end
    error('%d terms not supported\n', terms);
end

function [fun, thres] = smoothfun3_simple_taylor_setup(terms)
    fun = @(d, rd) smoothfun3_simple_taylor(d, rd, terms);
    if terms == 1
        thres = 2.98e-3;
    elseif terms == 2
        thres = 1.40e-2;
    elseif terms == 3
        thres = 3.69e-2;
    elseif terms == 4
        thres = 7.20e-2;
    else
        error('%d terms not supported\n', terms);
    end
end

% vim:set shiftwidth=4 softtabstop=4:
