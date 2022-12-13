function SE_warnings(key)
%SE_warnings  Turn on or off warnings for the SE package
%
%   SE_warnings('on')
%   SE_warnings('off')

warning(key, 'SE3P:SlowGridding');
warning(key, 'SE2P:SlowGridding');
warning(key, 'SE1P:SlowGridding');
warning(key, 'SE0P:SlowGridding');

warning(key, 'SE3P:PolynomialShapeFactor');
warning(key, 'SE2P:PolynomialShapeFactor');
warning(key, 'SE1P:PolynomialShapeFactor');
warning(key, 'SE0P:PolynomialShapeFactor');

warning(key, 'SE3P:ChangedGridRes');
warning(key, 'SE2P:ChangedGridRes');
warning(key, 'SE1P:ChangedGridRes');
