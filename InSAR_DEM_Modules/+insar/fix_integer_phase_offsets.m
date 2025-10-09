function [phzFixed, cycleMap] = fix_integer_phase_offsets(phzUnwrapped)
    % cycleMap = floor(phzUnwrapped / (2*pi));
    cycleMap = round((phzUnwrapped + pi) / (2*pi));
    phzFixed = phzUnwrapped - 2*pi * cycleMap;
end