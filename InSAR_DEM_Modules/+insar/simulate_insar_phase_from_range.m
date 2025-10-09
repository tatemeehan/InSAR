function [unwrappedPhz,wrappedPhz] = simulate_insar_phase_from_range(lambda, slantRange1, slantRange2)
% simulate insar phase from range
% Follows Equation 2.1 of InSAR Principles: Guidelines for SAR 
% Interferometry Processing and Interpretation (ESA TM-19) Part B
% https://www.esa.int/About_Us/ESA_Publications/InSAR_Principles_Guidelines_for_SAR_Interferometry_Processing_and_Interpretation_br_ESA_TM-19
% to Calcuate the Unwrapped interferrometric Phase of Bisatic InSAR

% Unwrapped Phase
unwrappedPhz = (4.*pi./lambda).*(slantRange2 - slantRange1);
% Wrapped Phase
wrappedPhz = mod(unwrappedPhz + pi,2.*pi) - pi;

end