% Input: Full y matrix y
%        Time-Wise index k
% Output: Which stations are active at any given time
function [stations] = FindStations(y,k)

  stations = (unique(ceil(find(abs(y(:,k)) > 0)/3)));

end