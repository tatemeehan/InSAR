function [I_corr, model] = remove_phzramp_geometry_gammaColumn(Ifilt, coh, geom, trajM, trajS, lambda, opts)
% Scale φ_FE per column using vertical (range) gradient WLS (wrap-safe).
% opts.cohThresh (0.30), opts.sigmaFE (30–80), opts.sigmaGamma (100–300), opts.useGX (false)

if nargin<7, opts = struct(); end
if ~isfield(opts,'cohThresh'),  opts.cohThresh  = 0.30; end
if ~isfield(opts,'sigmaFE'),    opts.sigmaFE    = 40;   end
if ~isfield(opts,'sigmaGamma'), opts.sigmaGamma = 150;  end
if ~isfield(opts,'useGX'),      opts.useGX      = false; end

k = 4*pi/lambda;
[M,N] = size(Ifilt);

% ----- geometry to get ŝ and B -----
X=geom.X; Y=geom.Y; Z=geom.Z;
idxM = double(geom.closestIndex_master);
idxS = double(geom.closestIndex_slave);
PM = trajM.pos;  PS = trajS.pos;  tM = (1:size(PM,1))'; tS = (1:size(PS,1))';

rmx = interp1(tM, PM(:,1), idxM, 'linear', NaN);
rmy = interp1(tM, PM(:,2), idxM, 'linear', NaN);
rmz = interp1(tM, PM(:,3), idxM, 'linear', NaN);
rsx = interp1(tS, PS(:,1), idxS, 'linear', NaN);
rsy = interp1(tS, PS(:,2), idxS, 'linear', NaN);
rsz = interp1(tS, PS(:,3), idxS, 'linear', NaN);

dx=X-rmx; dy=Y-rmy; dz=Z-rmz;  R=sqrt(dx.^2+dy.^2+dz.^2); R(R==0)=NaN;
sx=dx./R; sy=dy./R; sz=dz./R;                              % LOS
Bx=rsx-rmx; By=rsy-rmy; Bz=rsz-rmz;                        % baseline

phiFE = -k*(Bx.*sx + By.*sy + Bz.*sz);                     % unwrapped FE
phiFE(~isfinite(phiFE)) = NaN;
phiFE = nanGaussLP(phiFE, opts.sigmaFE);                   % gentle LP
phiFE = phiFE - median(phiFE(isfinite(phiFE)),'omitnan');  % gauge

% ----- wrap-safe grads from data & FE -----
[gx, gy, Wx, Wy] = phaseGradAndWeights(Ifilt, coh);        % chain-safe helper
[px, py]         = gradsFinite(phiFE);                     % FE grads

% ----- per-column γ(j) from vertical gradient (optionally add gx/px) -----
gamma = nan(1,N);
for j=1:N
    num=0; den=0;

    % vertical pairs at column j
    if j<=N
        mv = isfinite(gy(:,min(j,size(gy,2)))) & isfinite(Wy(:,min(j,size(Wy,2)))) & ...
             isfinite(py(:,min(j,size(py,2))));
        if any(mv)
            w  = Wy(mv, min(j,size(Wy,2)));
            yv = gy(mv, min(j,size(gy,2)));
            pv = py(mv, min(j,size(py,2)));
            num = num + nansum(w .* pv .* yv);
            den = den + nansum(w .* pv .* pv);
        end
    end

    % (optional) include horizontal pairs adjacent to j for stability
    if opts.useGX
        if j<=N-1
            mh = isfinite(gx(:,j)) & isfinite(Wx(:,j)) & isfinite(px(:,j));
            if any(mh)
                w  = Wx(mh,j);
                yh = gx(mh,j);
                ph = px(mh,j);
                num = num + nansum(w .* ph .* yh);
                den = num + nansum(w .* ph .* ph);
            end
        end
        if j>=2
            mh = isfinite(gx(:,j-1)) & isfinite(Wx(:,j-1)) & isfinite(px(:,j-1));
            if any(mh)
                w  = Wx(mh,j-1);
                yh = gx(mh,j-1);
                ph = px(mh,j-1);
                num = num + nansum(w .* ph .* yh);
                den = den + nansum(w .* ph .* ph);
            end
        end
    end

    if den>0, gamma(j)=num/den; end
end

% smooth γ along columns (NaN-aware)
gamma = smooth1D_nan(gamma(:), opts.sigmaGamma).';

% ----- predict and subtract -----
phi_pred = phiFE .* repmat(gamma, M, 1);
phi_pred = phi_pred - median(phi_pred(isfinite(phi_pred)),'omitnan');  % gauge
I_corr   = Ifilt .* exp(-1i*phi_pred);

% pack
model.gamma_j  = gamma;
model.phi_FE   = phiFE;
model.phi_pred = phi_pred;
end

% ---- helpers ----
function [gx, gy, Wx, Wy] = phaseGradAndWeights(I, coh)
% I   : MxN complex interferogram
% coh : MxN coherence
% gx  : Mx(N-1) horizontal wrapped gradient (rad)
% gy  : (M-1)xN vertical wrapped gradient (rad)
% Wx,Wy : pair weights = sqrt(coh_pair1 * coh_pair2)

    [M,N] = size(I);
    U  = I ./ max(abs(I), eps);

    % Preallocate
    gx = NaN(M, N-1);  gy = NaN(M-1, N);
    Wx = NaN(M, N-1);  Wy = NaN(M-1, N);

    % ----- Horizontal pairs (j vs j+1) -----
    U1x = U(:,1:N-1);           U2x = U(:,2:N);
    C1x = coh(:,1:N-1);         C2x = coh(:,2:N);
    mx  = isfinite(U1x) & isfinite(U2x) & isfinite(C1x) & isfinite(C2x);

    tmpGx = U2x .* conj(U1x);           % complex ratio
    tmpWx = sqrt(max(C1x,0) .* max(C2x,0));

    gx(mx) = angle(tmpGx(mx));
    Wx(mx) = tmpWx(mx);

    % ----- Vertical pairs (i vs i+1) -----
    U1y = U(1:M-1,:);           U2y = U(2:M,:);
    C1y = coh(1:M-1,:);         C2y = coh(2:M,:);
    my  = isfinite(U1y) & isfinite(U2y) & isfinite(C1y) & isfinite(C2y);

    tmpGy = U2y .* conj(U1y);
    tmpWy = sqrt(max(C1y,0) .* max(C2y,0));

    gy(my) = angle(tmpGy(my));
    Wy(my) = tmpWy(my);
end


function [px, py] = gradsFinite(A)
% A   : MxN real field (e.g., phi_FE after LP)
% px  : Mx(N-1) horizontal diff A(:,j+1)-A(:,j)
% py  : (M-1)xN vertical   diff A(i+1,:)-A(i,:)

    [M,N] = size(A);

    % Horizontal
    A1x = A(:,1:N-1);  A2x = A(:,2:N);
    mx  = isfinite(A1x) & isfinite(A2x);
    Dx  = A2x - A1x;
    px  = NaN(M, N-1);
    px(mx) = Dx(mx);

    % Vertical
    A1y = A(1:M-1,:);  A2y = A(2:M,:);
    my  = isfinite(A1y) & isfinite(A2y);
    Dy  = A2y - A1y;
    py  = NaN(M-1, N);
    py(my) = Dy(my);
end


function Aout = nanGaussLP(A, sigma)
if sigma<=0 || all(~isfinite(A(:))), Aout=A; return; end
L = max(3, ceil(5*sigma)); x = (-L:L);
g = exp(-(x.^2)/(2*sigma^2)); g = g/sum(g);
mask = isfinite(A); A(~mask)=0;
num = conv2(conv2(A, g, 'same'), g', 'same');
den = conv2(conv2(double(mask), g, 'same'), g', 'same');
Aout = num ./ max(den, eps); Aout(den==0)=NaN;
end

function b = smooth1D_nan(a, sigma)
if sigma<=0, b=a; return; end
L=max(3,ceil(5*sigma)); x=(-L:L); g=exp(-(x.^2)/(2*sigma^2)); g=g/sum(g);
m=isfinite(a); a2=a; a2(~m)=0; m2=double(m);
num=conv(a2,g,'same'); den=conv(m2,g,'same');
b=num./max(den,eps); b(den==0)=NaN;
end