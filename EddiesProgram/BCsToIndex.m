function M_indexes = BCsToIndex(M)
% remember this can be dodgy if it's only used (once) in initFilter, dont
% put this inside a loop because it will fuck efficiency

% Eddies order is:
%
% \
% |
% |
% |
% x

M_indexes =[];
for nn=1:size(M,2)
    for mm=1:size(M,1)
        if M(mm,nn) ~= 0
            M_indexes = [M_indexes, (nn-1)*size(M,1)+mm];
        end
    end
end