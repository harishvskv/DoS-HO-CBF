function T = time_to_floor(h, hdot, a_d, h_floor)
%TIME_TO_FLOOR  Predicted time for the barrier to reach the feasibility floor.

h    = h(:);
hdot = hdot(:);
if isscalar(h) && ~isscalar(hdot), h    = repmat(h,size(hdot)); end
if isscalar(hdot) && ~isscalar(h), hdot = repmat(hdot,size(h)); end

T = inf(size(h));

already = h <= h_floor;
T(already) = 0;

act = ~already;
if ~any(act), return; end

if abs(a_d) < eps
    lin = act & (hdot < 0);
    T(lin) = (h_floor - h(lin))./hdot(lin);
    return
end

disc = hdot.^2 + 2*a_d*(h - h_floor);
ok   = act & (disc >= 0);

Troot = (hdot(ok) + sqrt(disc(ok)))/a_d;
Troot(~isreal(Troot) | Troot < 0) = Inf;
T(ok) = Troot;
end
