function scene = build_insar_scene_from_ML(ML, LiDAR, opts)
%BUILD_INSAR_SCENE_FROM_ML Package ML snow/soil fields for InSAR modeling.

if nargin < 3, opts = struct(); end
if ~isfield(opts,'snowDepthField'), opts.snowDepthField = 'RF'; end
if ~isfield(opts,'soilLossTangent'), opts.soilLossTangent = 0.05; end
if ~isfield(opts,'f0'), opts.f0 = 1.3e9; end

c0 = 299792458;
scene.radar.f0 = opts.f0;
scene.radar.lambda0 = c0 / opts.f0;
scene.radar.k0 = 2*pi / scene.radar.lambda0;

% --- Snow dielectric field ---
scene.snow.eps = complex( ...
    ML.GPR.RealPermittivity.mean, ...
    ML.GPR.ImagPermittivity.mean);

scene.snow.epsReal = real(scene.snow.eps);
scene.snow.epsImag = imag(scene.snow.eps);

% --- Snow state fields ---
scene.snow.rhoWet = ML.GPR.Density.mean;
scene.snow.LWCpct = ML.GPR.LWC.mean;
scene.snow.LWC = scene.snow.LWCpct ./ 100;

if isfield(ML.GPR,'DensityDry') && isfield(ML.GPR.DensityDry,'mean')
    scene.snow.rhoDry = ML.GPR.DensityDry.mean;
else
    scene.snow.rhoDry = scene.snow.rhoWet - 10 .* scene.snow.LWCpct;
end

% --- Snow depth ---
switch upper(opts.snowDepthField)
    case 'RF'
        scene.snow.depth = LiDAR.A.RF;
    case 'RAW'
        scene.snow.depth = LiDAR.A.A;
    case 'GPR'
        scene.snow.depth = LiDAR.A.GPR;
    otherwise
        error('Unknown opts.snowDepthField: %s', opts.snowDepthField);
end

% --- Soil interface field ---
soilEr = ML.GPR.SoilPermittivity.mean;
soilEi = opts.soilLossTangent .* soilEr;

scene.soil.epsSurface = complex(soilEr, soilEi);
scene.soil.epsReal = soilEr;
scene.soil.epsImag = soilEi;
scene.soil.lossTangent = opts.soilLossTangent;

% --- Derived wave fields ---
scene.snow.n = sqrt(scene.snow.eps);
scene.snow.k = scene.radar.k0 .* scene.snow.n;
scene.snow.beta = real(scene.snow.k);
scene.snow.alpha = imag(scene.snow.k);

scene.soil.nSurface = sqrt(scene.soil.epsSurface);

% Normal-incidence first-pass reflection coefficient
scene.interface.Rsg = ...
    (scene.soil.nSurface - scene.snow.n) ./ ...
    (scene.soil.nSurface + scene.snow.n);

% First-pass propagation diagnostics
scene.phase.oneWaySnow = scene.snow.beta .* scene.snow.depth;
scene.phase.twoWaySnow = 2 .* scene.phase.oneWaySnow;

scene.amp.oneWayAtten = exp(-scene.snow.alpha .* scene.snow.depth);
scene.amp.twoWayAtten = exp(-2 .* scene.snow.alpha .* scene.snow.depth);

% Common valid mask
scene.mask.valid = isfinite(scene.snow.epsReal) & ...
                   isfinite(scene.snow.epsImag) & ...
                   isfinite(scene.snow.depth) & ...
                   isfinite(scene.soil.epsReal);

end