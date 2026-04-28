function dphi = phase_from_pass_pair(E1, E2)
% Return interferometric phase = angle(E1 .* conj(E2)) per channel.
E1 = ensure_2x2(E1);
E2 = ensure_2x2(E2);
dphi = angle(E1 .* conj(E2));  % elementwise, returns 2x2
end

function X2 = ensure_2x2(X)
if isempty(X), X2 = complex(zeros(2,2));
elseif isscalar(X), X2 = [X 0; 0 X];
else, X2 = X;
end
end
