function pairs = build_pairs(sarData, opts)
%BUILD_PAIRS Create list of interferometric pairs across and within dirs.
% pairs(p): struct('i1',i1,'b1',b1,'i2',i2,'b2',b2)
%
% opts.includeInter       (default true)
% opts.includeIntra       (default false)
% opts.burstMode          'sameIndex' (default) | 'allCombos'   % legacy
% opts.burstModeInter     'sameIndex' (default) | 'allCombos'
% opts.burstModeIntra     'neighbors' (default) | 'allCombos'
% opts.ensureIntraAllCombos (default true)  % augment neighbors -> all (b1<b2)
% opts.minBurstsPerDir    (default 1)
% opts.filterFn           @(i1,b1,i2,b2) -> true/false

if nargin < 2, opts = struct(); end
if ~isfield(opts,'includeInter'),         opts.includeInter = true; end
if ~isfield(opts,'includeIntra'),         opts.includeIntra = false; end
if ~isfield(opts,'burstMode'),            opts.burstMode = 'sameIndex'; end
if ~isfield(opts,'burstModeInter'),       opts.burstModeInter = opts.burstMode; end
if ~isfield(opts,'burstModeIntra'),       opts.burstModeIntra = 'neighbors'; end
if ~isfield(opts,'ensureIntraAllCombos'), opts.ensureIntraAllCombos = true; end
if ~isfield(opts,'minBurstsPerDir'),      opts.minBurstsPerDir = 1; end

pairs = struct('i1',{},'b1',{},'i2',{},'b2',{});
numDirs = numel(sarData);

    function try_add(i1,b1,i2,b2)
        if isfield(opts,'filterFn') && ~isempty(opts.filterFn)
            if ~opts.filterFn(i1,b1,i2,b2), return; end
        end
        pairs(end+1) = struct('i1',i1,'b1',b1,'i2',i2,'b2',b2); %#ok<AGROW>
    end

% -------- Inter-directory --------
if opts.includeInter && numDirs >= 2
    for i1 = 1:numDirs-1
        n1 = numel(sarData(i1).slc);
        if n1 < opts.minBurstsPerDir, continue; end
        for i2 = i1+1:numDirs
            n2 = numel(sarData(i2).slc);
            if n2 < opts.minBurstsPerDir, continue; end
            switch lower(opts.burstModeInter)
                case 'sameindex'
                    Nb = min(n1,n2);
                    for b = 1:Nb
                        try_add(i1,b,i2,b);
                    end
                case 'allcombos'
                    for b1 = 1:n1
                        for b2 = 1:n2
                            try_add(i1,b1,i2,b2);
                        end
                    end
                otherwise
                    error('Unknown burstModeInter: %s', opts.burstModeInter);
            end
        end
    end
end

% -------- Intra-directory --------
if opts.includeIntra
    for i1 = 1:numDirs
        n1 = numel(sarData(i1).slc);
        if n1 < 2 || n1 < opts.minBurstsPerDir, continue; end

        switch lower(opts.burstModeIntra)
            case {'neighbors','adjacent','sameindex'}
                for b1 = 1:n1-1
                    try_add(i1,b1,i1,b1+1);
                end
            case 'allcombos'
                for b1 = 1:n1-1
                    for b2 = b1+1:n1
                        try_add(i1,b1,i1,b2);
                    end
                end
            otherwise
                error('Unknown burstModeIntra: %s', opts.burstModeIntra);
        end

        % Augment to all (b1<b2) if requested
        if opts.ensureIntraAllCombos
            have = false(n1,n1);
            for p = 1:numel(pairs)
                if pairs(p).i1==i1 && pairs(p).i2==i1
                    have(pairs(p).b1, pairs(p).b2) = true;
                end
            end
            for b1 = 1:n1-1
                for b2 = b1+1:n1
                    if ~have(b1,b2)
                        try_add(i1,b1,i1,b2);
                    end
                end
            end
        end
    end
end
end

% function pairs = build_pairs(sarData, opts)
% %BUILD_PAIRS Create list of interferometric pairs across and within dirs.
% %   pairs(p) has fields: i1,b1,i2,b2
% %
% % opts.includeInter     (default true)   -> across directories
% % opts.includeIntra     (default false)  -> within each directory
% % opts.burstMode        'sameIndex' (default) | 'allCombos'   (legacy umbrella)
% % opts.burstModeInter   'sameIndex' (default) | 'allCombos'   (new)
% % opts.burstModeIntra   'neighbors' (default) | 'allCombos'   (new)
% % opts.ensureIntraAllCombos (default true) -> augment intra pairs to all (b1<b2)
% % opts.minBurstsPerDir  (default 1)
% % opts.filterFn         function handle @(i1,b1,i2,b2) -> true/false (optional)
% 
% if nargin < 2, opts = struct(); end
% if ~isfield(opts,'includeInter'),        opts.includeInter = true; end
% if ~isfield(opts,'includeIntra'),        opts.includeIntra = true; end
% if ~isfield(opts,'burstMode'),           opts.burstMode = 'sameIndex'; end
% if ~isfield(opts,'burstModeInter'),      opts.burstModeInter = opts.burstMode; end
% if ~isfield(opts,'burstModeIntra'),      opts.burstModeIntra = 'neighbors'; end
% if ~isfield(opts,'ensureIntraAllCombos'),opts.ensureIntraAllCombos = true; end
% if ~isfield(opts,'minBurstsPerDir'),     opts.minBurstsPerDir = 1; end
% 
% pairs = struct('i1',{},'b1',{},'i2',{},'b2',{});
% numDirs = numel(sarData);
% 
%     function try_add(i1,b1,i2,b2)
%         if isfield(opts,'filterFn') && ~isempty(opts.filterFn)
%             if ~opts.filterFn(i1,b1,i2,b2), return; end
%         end
%         pairs(end+1) = struct('i1',i1,'b1',b1,'i2',i2,'b2',b2); %#ok<AGROW>
%     end
% 
% %% Inter-directory pairs
% if opts.includeInter && numDirs >= 2
%     for i1 = 1:numDirs-1
%         n1 = numel(sarData(i1).slc);
%         if n1 < opts.minBurstsPerDir, continue; end
%         for i2 = i1+1:numDirs
%             n2 = numel(sarData(i2).slc);
%             if n2 < opts.minBurstsPerDir, continue; end
% 
%             switch lower(opts.burstModeInter)
%                 case 'sameindex'
%                     Nb = min(n1,n2);
%                     for b = 1:Nb
%                         try_add(i1,b,i2,b);
%                     end
%                 case 'allcombos'
%                     for b1 = 1:n1
%                         for b2 = 1:n2
%                             try_add(i1,b1,i2,b2);
%                         end
%                     end
%                 otherwise
%                     error('Unknown burstModeInter: %s', opts.burstModeInter);
%             end
%         end
%     end
% end
% 
% %% Intra-directory pairs
% if opts.includeIntra
%     for i1 = 1:numDirs
%         n1 = numel(sarData(i1).slc);
%         if n1 < 2 || n1 < opts.minBurstsPerDir, continue; end
% 
%         switch lower(opts.burstModeIntra)
%             case {'neighbors','adjacent','sameindex'}  % treat 'sameIndex' as neighbors intra-dir
%                 for b1 = 1:n1-1
%                     try_add(i1,b1,i1,b1+1);
%                 end
%             case 'allcombos'
%                 for b1 = 1:n1-1
%                     for b2 = b1+1:n1
%                         try_add(i1,b1,i1,b2);
%                     end
%                 end
%             otherwise
%                 error('Unknown burstModeIntra: %s', opts.burstModeIntra);
%         end
% 
%         % --- Augment to all combos if requested (guarantees (1,3), etc.) ---
%         if opts.ensureIntraAllCombos
%             % collect existing intra pairs for this dir
%             have = false(n1,n1);
%             for p = 1:numel(pairs)
%                 if pairs(p).i1==i1 && pairs(p).i2==i1
%                     have(pairs(p).b1, pairs(p).b2) = true;
%                 end
%             end
%             % add any missing (b1<b2)
%             for b1 = 1:n1-1
%                 for b2 = b1+1:n1
%                     if ~have(b1,b2)
%                         try_add(i1,b1,i1,b2);
%                     end
%                 end
%             end
%         end
%     end
% end
% end
% 
% % function pairs = build_pairs(sarData, opts)
% % %BUILD_PAIRS Create list of interferometric pairs across and within dirs.
% % % pairs(p) has fields: i1,b1,i2,b2
% % %
% % % opts.includeInter  (default true)  -> across directories
% % % opts.includeIntra  (default false) -> within each directory
% % % opts.burstMode     'sameIndex' (default) | 'allCombos'
% % % opts.minBurstsPerDir (default 1)
% % % opts.filterFn      function handle @(i1,b1,i2,b2) -> true/false (optional)
% % 
% % if nargin < 2, opts = struct(); end
% % if ~isfield(opts,'includeInter'), opts.includeInter = true; end
% % if ~isfield(opts,'includeIntra'), opts.includeIntra = false; end
% % if ~isfield(opts,'burstMode'),    opts.burstMode    = 'sameIndex'; end
% % if ~isfield(opts,'minBurstsPerDir'), opts.minBurstsPerDir = 1; end
% % 
% % pairs = struct('i1',{},'b1',{},'i2',{},'b2',{});
% % numDirs = numel(sarData);
% % 
% % % helper to add pair with optional user filter
% %     function try_add(i1,b1,i2,b2)
% %         if isfield(opts,'filterFn') && ~isempty(opts.filterFn)
% %             if ~opts.filterFn(i1,b1,i2,b2), return; end
% %         end
% %         pairs(end+1) = struct('i1',i1,'b1',b1,'i2',i2,'b2',b2); %#ok<AGROW>
% %     end
% % 
% % % Inter-directory pairs
% % if opts.includeInter && numDirs >= 2
% %     for i1 = 1:numDirs-1
% %         n1 = numel(sarData(i1).slc);
% %         if n1 < opts.minBurstsPerDir, continue; end
% %         for i2 = i1+1:numDirs
% %             n2 = numel(sarData(i2).slc);
% %             if n2 < opts.minBurstsPerDir, continue; end
% % 
% %             switch lower(opts.burstMode)
% %                 case 'sameindex'
% %                     Nb = min(n1,n2);
% %                     for b = 1:Nb
% %                         try_add(i1,b,i2,b);
% %                     end
% %                 case 'allcombos'
% %                     for b1 = 1:n1
% %                         for b2 = 1:n2
% %                             try_add(i1,b1,i2,b2);
% %                         end
% %                     end
% %                 otherwise
% %                     error('Unknown burstMode: %s', opts.burstMode);
% %             end
% %         end
% %     end
% % end
% % 
% % % Intra-directory pairs
% % if opts.includeIntra
% %     for i1 = 1:numDirs
% %         n1 = numel(sarData(i1).slc);
% %         if n1 < 2 || n1 < opts.minBurstsPerDir, continue; end
% % 
% %         switch lower(opts.burstMode)
% %             case 'sameindex'
% %                 % sameIndex doesn’t make sense intra-dir unless you want (b,b); skip
% %                 % Instead, pair neighbors by index:
% %                 for b1 = 1:n1-1
% %                     b2 = b1+1;
% %                     try_add(i1,b1,i1,b2);
% %                 end
% %             case 'allcombos'
% %                 for b1 = 1:n1-1
% %                     for b2 = b1+1:n1
% %                         try_add(i1,b1,i1,b2);
% %                     end
% %                 end
% %             otherwise
% %                 error('Unknown burstMode: %s', opts.burstMode);
% %         end
% %     end
% % end
% % end
