function out = combineWeightEntries(e1,e2)
% Union indices, average values (normalize). Empty-safe.
if isempty(e1) && isempty(e2), out = []; return; end
if isempty(e1), out = e2; return; end
if isempty(e2), out = e1; return; end
idx = unique([e1.indices(:); e2.indices(:)]);
val = zeros(size(idx));
[~,l1] = ismember(e1.indices(:), idx); val(l1) = val(l1) + e1.value(:);
[~,l2] = ismember(e2.indices(:), idx); val(l2) = val(l2) + e2.value(:);
s = sum(val); if s>0, val = val/s; end
out.indices = idx; out.value = val;
end
