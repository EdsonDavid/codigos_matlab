function r = radio(x)
if x < 1
    r = -3*x + 4.5;
else
    r = repmat(1.5,size(x));
end
