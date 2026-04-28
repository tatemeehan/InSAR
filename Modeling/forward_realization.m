function out = forward_realization(theta_i,f,Hs,eps_snow,eps_g0,pars)
% forward_realization: assumes phases are already set by caller (e.g., coherence_mc).
% If you want random phases, set them in the caller or add a separate helper.

% out = insar_forward_snow_soil_benchmark(theta_i,f,Hs,eps_snow,eps_g0,pars);
% Call 4 component Model
out = insar_forward_snow_soil_4component(theta_i, f, Hs, eps_snow, eps_g0, pars);
end

% function out = forward_realization(theta_i,f,Hs,eps_snow,eps_g0,pars)
% % Randomize unresolved phases only if caller hasn't set them.
% 
% if isfield(pars,'rough') && ~isempty(pars.rough)
%     if isfield(pars.rough,'Cdiff') && pars.rough.Cdiff > 0
%         if ~isfield(pars.rough,'psi') || isempty(pars.rough.psi)
%             pars.rough.psi = 2*pi*rand;
%         end
%     end
%     if isfield(pars.rough,'Xpol') && pars.rough.Xpol > 0
%         if ~isfield(pars.rough,'psi_x') || isempty(pars.rough.psi_x)
%             pars.rough.psi_x = 2*pi*rand;
%         end
%     end
% end
% 
% if isfield(pars,'snow') && ~isempty(pars.snow)
%     if isfield(pars.snow,'Xpol') && pars.snow.Xpol > 0
%         if ~isfield(pars.snow,'psi_x') || isempty(pars.snow.psi_x)
%             pars.snow.psi_x = 2*pi*rand;
%         end
%     end
% end
% 
% out = insar_forward_snow_soil_benchmark(theta_i,f,Hs,eps_snow,eps_g0,pars);
% end

% function out = forward_realization(theta_i,f,Hs,eps_snow,eps_g0,pars)
% % Randomize unresolved phases (one realization)
% if isfield(pars,'rough') && isfield(pars.rough,'Cdiff') && pars.rough.Cdiff > 0
%     if ~isfield(pars.rough,'psi') || isempty(pars.rough.psi)
%         pars.rough.psi = 2*pi*rand;
%     end
% end
% 
% if isfield(pars,'rough') && isfield(pars.rough,'Xpol') && pars.rough.Xpol > 0
%     if ~isfield(pars.rough,'psi_x') || isempty(pars.rough.psi_x)
%         pars.rough.psi_x = 2*pi*rand;
%     end
% end
% 
% if isfield(pars,'snow') && isfield(pars.snow,'Xpol') && pars.snow.Xpol > 0
%     if ~isfield(pars.snow,'psi_x') || isempty(pars.snow.psi_x)
%         pars.snow.psi_x = 2*pi*rand;
%     end
% end
% 
% out = insar_forward_snow_soil_benchmark(theta_i,f,Hs,eps_snow,eps_g0,pars);
% end
