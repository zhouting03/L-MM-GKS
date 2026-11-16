function root = setup_lplq()
% 添加本代码包的路径，不修改用户已有 MMGKS/LMMGKS 文件。
root = fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'solver'),fullfile(root,'operators'), ...
    fullfile(root,'experiments'),fullfile(root,'tests'));
end
