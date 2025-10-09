function S_slave_rephased = rephase_slave_to_master(S_slave, deltaR, lambda)
% S_slave:  complex slave SLC coregistered to master grid
% deltaR :  r2 - r1 map (meters) from your trajectory → pixel geometry
% lambda :  wavelength (meters), e.g., L-band ≈ 0.238 m
%
% Output is a *new* slave SLC whose phase has master range history.

% Most SAR processors store phase as exp(-j*4π r / λ).
% To change r2 → r1 multiply by exp(-j*4π (r1 - r2)/λ) = exp(+j*4π (r2 - r1)/λ)
S_slave_rephased = S_slave .* exp( 1i * (4*pi/lambda) * deltaR );
end
