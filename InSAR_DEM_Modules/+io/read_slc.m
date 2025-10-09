function slc = read_slc(filename, frameSize, nanMask)
%READ_SLC Load complex SLC binary data, swap bytes, and apply NaN mask
%
%   slc = read_slc(filename, frameSize, nanMask)
%
%   Inputs:
%       filename  - full path to the .slc file
%       frameSize - [rows, cols] of the image
%       nanMask   - logical mask for invalid regions
%
%   Output:
%       slc - complex SLC matrix with NaNs applied

    reader = dsp.BinaryFileReader(filename, ...
        'SamplesPerFrame', frameSize(1), ...
        'NumChannels', frameSize(2), ...
        'DataType', 'single', ...
        'IsDataComplex', true);

    slc = swapbytes(reader());
    if nargin > 2 && ~isempty(nanMask)
        slc(nanMask) = NaN;
    end
end