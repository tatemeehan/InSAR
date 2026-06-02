function stats = pearson_phase_components(phiData, phiModel, mask)

valid = mask & isfinite(phiData) & isfinite(phiModel);

x = double(phiData(valid));
y = double(phiModel(valid));

cu = cos(x);
su = sin(x);

cm = cos(y);
sm = sin(y);

Rc = corrcoef(cu, cm, 'Rows', 'complete');
Rs = corrcoef(su, sm, 'Rows', 'complete');

stats.r_cos = Rc(1,2);
stats.r_sin = Rs(1,2);
stats.r_mean = mean([stats.r_cos, stats.r_sin], 'omitnan');

end