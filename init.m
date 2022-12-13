function init()

root = fileparts(mfilename('fullpath')); % path to the directory containing this file
setenv('SE_reset_ROOT', root);

addpath([root, '/bin']);
addpath([root, '/SE_common']);
addpath([root, '/SE3P/src']);
addpath([root, '/SE2P/src']);
addpath([root, '/SE1P/src']);
addpath([root, '/SE1P/src/struve']);
addpath([root, '/SE0P/src']);

end
