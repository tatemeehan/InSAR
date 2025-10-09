function out = nanmedfilt2(in, windowSize)
%NANMEDFILT2 Median filter that ignores NaNs
%
%   out = nanmedfilt2(in, [m n]) applies a median filter of size [m n]
%   while ignoring NaN values in the neighborhood.
%
%   Example: out = nanmedfilt2(data, [3 3]);

if numel(windowSize) == 1
    windowSize = [windowSize windowSize];
end

padSize = floor(windowSize / 2);
inPad = padarray(in, padSize, NaN, 'both');
out = NaN(size(in));

for i = 1:size(in,1)
    for j = 1:size(in,2)
        window = inPad(i:i+windowSize(1)-1, j:j+windowSize(2)-1);
        validVals = window(~isnan(window));
        if ~isempty(validVals)
            out(i,j) = median(validVals);
        end
    end
end
end
