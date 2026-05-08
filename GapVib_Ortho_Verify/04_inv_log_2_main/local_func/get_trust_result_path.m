function trustPath = get_trust_result_path(thisDir, stem)
%get_trust_result_path  Resolve a saved trust-domain result path.

if nargin < 2 || isempty(stem)
    stem = "inv_log_2_trust_domain";
end
trustDir = fullfile(thisDir, 'trust_domain_result');
if ~exist(trustDir, 'dir')
    mkdir(trustDir);
end
trustPath = fullfile(trustDir, char(stem + ".mat"));
end
