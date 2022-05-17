FUNCTION error_ratio, x, y, xerr,yerr
return, sqrt((xerr/y)^2 + (x*yerr/y^2)^2)
END
