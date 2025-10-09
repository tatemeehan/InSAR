function tf = isSameGeometry(r1, r2, r1prev, r2prev, tol)
%ISSAMEGEOMETRY Check if (r1, r2) is geometrically equivalent to (r1prev, r2prev)
d1 = mean(vecnorm(r1 - r1prev, 2, 2));
d2 = mean(vecnorm(r2 - r2prev, 2, 2));
d_rev = mean(vecnorm(r1 - r2prev, 2, 2)) + mean(vecnorm(r2 - r1prev, 2, 2));
tf = (d1 < tol && d2 < tol) || (d_rev < 2*tol);
end