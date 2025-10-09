function [I_corr, phi_pred_back, thetaDeg, mdl] = ...
    remove_phzramp_poisson_rot(Ifilt, coh, headingDeg, sigmaRange, sigmaAz, lambda, cohThresh)

if nargin<4||isempty(sigmaRange), sigmaRange = 100;  end
if nargin<5||isempty(sigmaAz),    sigmaAz    = 25; end
if nargin<6||isempty(lambda),     lambda     = 1e-3; end
if nargin<7||isempty(cohThresh),  cohThresh  = 0.30; end

thetaDeg = headingDeg;   % see “Angle sanity check” below

% rotate into azimuth–range frame
If_rot  = rotfill(Ifilt, -thetaDeg, 'bilinear', 'crop', NaN);
coh_rot = rotfill(coh,   -thetaDeg, 'bilinear', 'crop', NaN);

% anisotropic Poisson in rotated frame (columns≈range, rows≈azimuth)
% [Icorr_rot, mdl] = insar.remove_phzramp_poisson_aniso( ...
%     If_rot, coh_rot, sigmaRange, sigmaAz, lambda, [], 0, 2, cohThresh);
[Icorr_rot, mdl] = insar.remove_phzramp_multiscale(If_rot, coh_rot);

% rotate the predicted screen back and subtract on the original grid
phi_pred_rot  = mdl;%angle(Icorr_rot);%mdl.phi_pred;
phi_pred_back = rotfill(phi_pred_rot, +thetaDeg, 'bilinear', 'crop', NaN);
I_corr        = Ifilt .* exp(-1i * phi_pred_back);
end



function Arot = rotfill(A, angDeg, method, bbox, fillVal)
% Rotate array A by angDeg (CCW), keep size ('crop'), and fill newly
% exposed pixels with fillVal (NaN). Works for real or complex.
if nargin<3||isempty(method), method='bilinear'; end
if nargin<4||isempty(bbox),   bbox='crop';      end
if nargin<5, fillVal = NaN; end

% rotate data
if isreal(A)
    Arot = imrotate(A, angDeg, method, bbox);
else
    Arot = imrotate(real(A), angDeg, method, bbox) ...
         + 1i*imrotate(imag(A), angDeg, method, bbox);
end

% build a mask of newly created pixels and set them to fillVal
M = imrotate(ones(size(A),'single'), angDeg, method, bbox);
edge = (M < 0.99);     % robust threshold
if isreal(Arot)
    Arot(edge) = fillVal;
else
    Arot(edge) = complex(fillVal, fillVal);
end
end
