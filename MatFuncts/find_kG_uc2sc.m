function [k_sc, G_sc] = find_kG_uc2sc(k_uc,M_uc2sc)
% This script is to convert the k point coordinate (direct corrdinate) in
% the primitive cell (large BZ) into k point corrdinate (direct corrdinate)
% in the sueprcell corresponded BZ (small BZ), and link them with a G
% vector that is defined for supercell
% Here k_sc and k_uc are in (1x3) row vector
% M_uc2sc is a (3x3) matrix that makes the transformation from unit cell to supercell
   k_sc = M_uc2sc * k_uc';
   k_sc = k_sc';
   G_sc = round(k_sc);
   k_sc = k_sc - G_sc; % so that the k_sc is in the range [-0.5, 0.5]
end
