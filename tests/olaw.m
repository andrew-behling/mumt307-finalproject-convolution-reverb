function olaw( type, N, olap, nWindows )
% OLAW  Plot overlapping windows
%
%   OLAW( TYPE, N, OLAP, nWINDOWS) plots nWINDOWS, each of length N,
%   overlapped by the percentage OLAP.  The window type is
%   specified with the string variable TYPE, which can be either
%   'triangle', 'hann', 'hamm', or 'blackman'. Note that there may be
%   problems in achieving constant overlap due to even/odd values of N and
%   the values of the windows at boundaries.
%
%   By Gary P. Scavone, McGill University, 2007.

if ( nargin ~= 4 )
  error('Number of arguments is incorrect.');
  return
end

hop = floor( N * (100-olap) / 100 );
nY = N + (nWindows-1) * hop;

if strcmp( type, 'triangle' )
  w = window( @triang, N );
elseif strcmp( type, 'hann' )
  w = window( @hann, N );
elseif strcmp( type, 'hamm' )
  w = window( @hamming, N );
elseif strcmp( type, 'blackman' )
  w = window( @blackman, N );
elseif strcmp( type, 'flattop')
  w = window( @ flattopwin, N);
else
  error('Window type argument is incorrect.');
  return;
end

y = zeros( nY, 1 );
y(1:N) = w;
for n = 2 : nWindows
  iStart = (n-1) * hop;
  y( iStart+1:iStart+N) = y( iStart+1:iStart+N) + w;
end

plot( y );
xlim( [1 nY] )


