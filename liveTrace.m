function liveTrace(p)
% liveTrace - adds the TCP position p = [x y z] (m) to the live figure of projectA_main
% Called during the Simscape simulation by a MATLAB Function block:
%     function liveTraceBlock(p)
%     coder.extrinsic('liveTrace');
%     liveTrace(p);
%     end
persistent h3 h2
if isempty(h3) || ~isvalid(h3(1)) || isempty(h2) || ~isvalid(h2(1))
    h3 = findobj(0, 'Tag', 'TCPtrace3D');                 % created in projectA_main, section 10
    h2 = findobj(0, 'Tag', 'TCPtrace2D');
    if isempty(h3) || isempty(h2), return; end
end
p = double(p(:));
addpoints(h3(1), p(1), p(2), p(3));                       % full path: approach, letter, return
ud = h2(1).UserData;
if abs(ud.rS(p(3)) - hypot(p(1), p(2))) < ud.tol          % pen on the wall -> draw on the unrolled view
    addpoints(h2(1), ud.R0*atan2(p(2), p(1)), p(3));
end
drawnow limitrate;
end