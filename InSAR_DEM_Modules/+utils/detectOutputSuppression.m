function isSuppressed = detectOutputSuppression(nout)
%DETECTOUTPUTSUPPRESSION Detect which outputs are suppressed with ~
%   isSuppressed = detectOutputSuppression(nout)
%
%   Returns a logical vector [true false ...] where true indicates that
%   output was suppressed using ~. If detection fails, it returns false(nout,1)

    try
        % Look at caller
        stack = dbstack('-completenames');
        if length(stack) < 2
            isSuppressed = false(1, nout);
            return;
        end

        callerFile = stack(2).file;
        callerLine = stack(2).line;

        % Read source code
        fileText = fileread(callerFile);
        lines = splitlines(fileText);
        callLine = strtrim(lines{callerLine});

        % Look for output pattern on this line
        tokens = regexp(callLine, '^\[([^\]]+)\]\s*=\s*', 'tokens');
        if isempty(tokens)
            isSuppressed = false(1, nout);
            return;
        end

        outList = strtrim(strsplit(tokens{1}{1}, ','));
        isSuppressed = cellfun(@(x) strcmp(x, '~'), outList);
        
        % Pad with false if fewer outputs than expected
        if length(isSuppressed) < nout
            isSuppressed(end+1:nout) = false;
        end

    catch
        % Fallback if anything breaks
        isSuppressed = false(1, nout);
    end
end
