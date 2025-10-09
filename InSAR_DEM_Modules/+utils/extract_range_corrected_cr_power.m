function CR_power = extract_range_corrected_cr_power(slc, slantRange, crIdxs, weights)
% EXTRACT_RANGE_CORRECTED_CR_POWER Extract power with range correction
amp = abs(slc);
amp_cr = amp(crIdxs);
range_cr = slantRange(crIdxs);

% % Range correction (linear power = A^2 * R^4)
% linPower = amp_cr.^2 .* (range_cr.^4);
% Range correction (linear power = A^2 * R^2)
linPower = amp_cr.^2 .* (range_cr.^2);
% No Range Correction Applied!!!
% linPower = amp_cr.^2;
CR_power = sum(linPower .* weights);
end