#include "..\Header Files\ssamporientimg.hpp"

/**
Matlab equivalent code:
% ssamporientimg - Realiza sub-sampling
%
% Usage:   ssamporientimg(orient, spacing)
%
%        orientim - Ridge orientation image (obtained from RIDGEORIENT)
%        spacing  - Sub-sampling interval to be used in ploting the
%                   orientation data the (Plotting every point is


function s_orient = ssamporientimg(orient, spacing)

if fix(spacing) ~= spacing
error('spacing must be an integer');
end

[rows, cols] = size(orient);

lw = 2;             % linewidth
len = 0.8*spacing;  % length of orientation lines

% Subsample the orientation data according to the specified spacing
s_orient = orient(spacing:spacing:rows-spacing, ...
spacing:spacing:cols-spacing);
return

*/