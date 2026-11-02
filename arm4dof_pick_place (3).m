%% ======================================================================
%  4-DOF R-R-R-R ARM FOR ROCKET PAYLOAD MODULE ASSEMBLY  (GNU Octave)
%  Picks 4 components from a rack on the arm's RIGHT and places each one in
%  its own numbered slot (S1-S4) in the payload bay on its LEFT.
%
%  Covers: 1 specs + 3D model   2 DH kinematics (FK/IK + validation)
%          3 path planning      4 motor + dynamics (fixed load)
%          5 ON-OFF / P / PI / PID control comparison
%
%  Run:  arm4dof_pick_place          (Octave 6 or newer, no packages needed)
%
%  Frame: base at origin, z up, arm "faces" +x.
%         Looking along +x, RIGHT = -y, LEFT = +y.
%  Units: SI inside the code (m, rad, kg, N*m); printed tables use mm / deg.
%% ======================================================================
1;   % marks this file as a script (functions follow)

% ---------------------------------------------------------------- params
function P = params()
  % DH table of the presentation: b1 (base height), a1 = 0, a2 = 300 mm, a3 = 250 mm, a4 = 0.
  % The end effector (frame 5) is mounted directly on the wrist joint: zero link length.
  P.a1 = 0;  P.a2 = 0.300;  P.a3 = 0.250;  P.a4 = 0;
  P.b1 = 0.150;                         % ASSUMED base height b1 (value not given in the slides)
  % Mass model (ASSUMED values - replace with your CAD/datasheet numbers)
  P.m2 = 1.2;  P.m3 = 0.9;              % link 2, link 3 (kg)
  P.me = 0.8;                           % elbow motor + gearbox (kg)
  P.mw = 0.4;                           % wrist motor + gripper body, placed at the wrist (kg)
  P.mL = 1.0;                           % FIXED payload / component (kg)
  P.I2 = P.m2*P.a2^2/12;  P.I3 = P.m3*P.a3^2/12;   % rod inertias
  P.Iw = 0.003;                         % wrist + gripper + payload inertia about the pitch axis (kg m^2)
  P.I1 = 0.01;                          % base column inertia about yaw axis
  P.g  = 9.81;
  % Joint limits
  P.qmin = deg2rad([-150; -30; -150; -150]);
  P.qmax = deg2rad([ 150; 135;  -10;  150]);
end

function Mo = motor_params()
  % One DC motor + gearbox per joint (same motor, different ratio)
  Mo.Ra = 1.0;  Mo.La = 0.5e-3;         % ohm, H
  Mo.Kt = 0.05; Mo.Ke = 0.05;           % N*m/A, V*s/rad
  Mo.Jm = 2e-5; Mo.Bm = 1e-5;           % rotor inertia, viscous friction
  Mo.Nr = [100; 120; 100; 50];          % gear ratios
  Mo.eta = 0.8;                         % gearbox efficiency
  Mo.Vmax = 24; Mo.Imax = 10;           % supply / current limit
  Mo.Von = 12;  Mo.db = deg2rad(0.5);   % ON-OFF drive voltage and dead-band
end

% ------------------------------------------------------------ kinematics
function T = dh(th, b, a, al)
  % Standard DH: Rz(theta) Tz(b) Tx(a) Rx(alpha)
  ct = cos(th); st = sin(th); ca = cos(al); sa = sin(al);
  T = [ct, -st*ca,  st*sa, a*ct;
       st,  ct*ca, -ct*sa, a*st;
        0,     sa,     ca,    b;
        0,      0,      0,    1];
end

function T = fk_dh(q, P)
  T12 = dh(q(1), P.b1, P.a1, pi/2);     % frame 1 -> 2  (base yaw)
  T23 = dh(q(2), 0,    P.a2, 0);        % frame 2 -> 3  (shoulder)
  T34 = dh(q(3), 0,    P.a3, 0);        % frame 3 -> 4  (elbow)
  T45 = dh(q(4), 0,    P.a4, 0);        % frame 4 -> 5  (wrist, end effector)
  T = T12*T23*T34*T45;                  % T1-5
end

function p = fk_closed(q, P)           % vectorised, q is 4 x n
  q23 = q(2,:) + q(3,:);  q234 = q23 + q(4,:);
  rho = P.a1 + P.a2*cos(q(2,:)) + P.a3*cos(q23) + P.a4*cos(q234);
  p = [rho.*cos(q(1,:));
       rho.*sin(q(1,:));
       P.b1 + P.a2*sin(q(2,:)) + P.a3*sin(q23) + P.a4*sin(q234)];
end

function q = ik(x, y, z, phi, P)       % vectorised, elbow-up; (x,y,z) = end effector, phi = tool pitch
  th1 = atan2(y, x);
  R = hypot(x, y) - P.a4*cos(phi) - P.a1;   % step back along link 4 to the wrist point
  s = z - P.b1 - P.a4*sin(phi);
  c3 = (R.^2 + s.^2 - P.a2^2 - P.a3^2) / (2*P.a2*P.a3);
  bad = abs(c3) > 1;
  c3(bad) = 1;
  th3 = -acos(c3);                     % negative = elbow up
  th2 = atan2(s, R) - atan2(P.a3*sin(th3), P.a2 + P.a3*cos(th3));
  th4 = phi - th2 - th3;
  q = [th1; th2; th3; th4];
  q(:, bad) = NaN;
end

function ok = within_limits(q, P)
  ok = all(q >= P.qmin - 1e-9 & q <= P.qmax + 1e-9, 1) & ~any(isnan(q), 1);
end

% ------------------------------------------------------------ trajectory
function s = trap(tau, f)              % normalised CUBIC polynomial profile (rest-to-rest)
  % s(tau) = 3 tau^2 - 2 tau^3 :  s(0)=0, s(1)=1, s'(0)=s'(1)=0.
  % (name and argument f kept from the old trapezoidal version so make_traj needs no change; f is unused)
  s = 3*tau.^2 - 2*tau.^3;
end

function Tr = make_traj(W, segs, phiW, P, dt, f)
  % segs rows: [type from to duration]
  %   type 1 = joint-space PTP, 2 = Cartesian straight line, 0 = dwell
  qW = ik(W(1,:), W(2,:), W(3,:), phiW, P);
  q = qW(:, segs(1,2));
  iend = zeros(1, size(segs,1));
  for k = 1:size(segs,1)
    ty = segs(k,1); a = segs(k,2); b = segs(k,3); T = segs(k,4);
    n = round(T/dt);  tau = (1:n)/n;  sp = trap(tau, f);
    if ty == 2
      pos = W(:,a) + (W(:,b)-W(:,a))*sp;
      qs = ik(pos(1,:), pos(2,:), pos(3,:), phiW(a), P);
    else
      qs = qW(:,a) + (qW(:,b)-qW(:,a))*sp;
    end
    q = [q, qs];
    iend(k) = size(q,2);
  end
  N = size(q,2);
  qd = zeros(4,N); qdd = zeros(4,N);
  for i = 1:4
    qd(i,:)  = gradient(q(i,:), dt);
    qdd(i,:) = gradient(qd(i,:), dt);
  end
  Tr.t = (0:N-1)*dt;  Tr.q = q;  Tr.qd = qd;  Tr.qdd = qdd;
  Tr.iend = iend;  Tr.qW = qW;
end

% -------------------------------------------------------------- dynamics
function [M, G] = MG(q, P)
  % Lumped masses: link-2 COM, link-3 COM, elbow motor, link-4 COM, payload at the TCP
  K  = [0.5 0 0; 1 0.5 0; 1 0 0; 1 1 0.5; 1 1 1];   % fractions along links 2, 3, 4
  mm = [P.m2; P.m3; P.me; P.mw; P.mL];
  q2 = q(2); q23 = q(2) + q(3); q234 = q23 + q(4);
  M = zeros(4);  G = zeros(4,1);
  for k = 1:5
    k2 = K(k,1); k3 = K(k,2); k4 = K(k,3);
    rho = P.a1 + k2*P.a2*cos(q2) + k3*P.a3*cos(q23) + k4*P.a4*cos(q234);
    dr = [-k2*P.a2*sin(q2) - k3*P.a3*sin(q23) - k4*P.a4*sin(q234);
          -k3*P.a3*sin(q23) - k4*P.a4*sin(q234);
          -k4*P.a4*sin(q234)];
    dz = [ k2*P.a2*cos(q2) + k3*P.a3*cos(q23) + k4*P.a4*cos(q234);
           k3*P.a3*cos(q23) + k4*P.a4*cos(q234);
           k4*P.a4*cos(q234)];
    M(1,1) = M(1,1) + mm(k)*rho^2;
    M(2:4,2:4) = M(2:4,2:4) + mm(k)*(dr*dr' + dz*dz');
    G(2:4) = G(2:4) + mm(k)*P.g*dz;
  end
  ev = [0 1 0 0]; M = M + P.I2*(ev'*ev);
  ev = [0 1 1 0]; M = M + P.I3*(ev'*ev);
  ev = [0 1 1 1]; M = M + P.Iw*(ev'*ev);
  M(1,1) = M(1,1) + P.I1;
end

function [tau, M] = invdyn(q, qd, qdd, P)
  [M, G] = MG(q, P);
  h = 1e-6;  dM = zeros(4,4,4);
  for k = 2:4                          % M depends on q2, q3, q4
    dq = zeros(4,1); dq(k) = h;
    [Mp, ~] = MG(q+dq, P);  [Mm, ~] = MG(q-dq, P);
    dM(:,:,k) = (Mp - Mm)/(2*h);
  end
  c = zeros(4,1);
  for i = 1:4                          % Christoffel symbols
    A = reshape(dM(i,:,:), 4, 4);
    B = dM(:,:,i);
    c(i) = qd' * (0.5*(A + A' - B)) * qd;
  end
  tau = M*qdd + c + G;
end

% --------------------------------------------------------- control sim
function out = run_ctrl(ctrl, S)
  Mo = S.Mo;  dt = S.dt;  N = numel(S.t);  ns = S.ns;  Ts = ns*dt;
  Nr = Mo.Nr;
  q = zeros(4,N);  Ilog = zeros(4,N);  Vlog = zeros(4,N);
  q(:,1) = S.qref(:,1);
  qm = Nr.*S.qref(:,1);  wm = zeros(4,1);  cur = zeros(4,1);
  integ = zeros(4,1);  ep = zeros(4,1);  u = zeros(4,1);
  for k = 1:N-1
    if mod(k-1, ns) == 0               % digital controller, Ts = 1 ms
      e = S.qref(:,k) - qm./Nr;
      switch ctrl
        case 'ONOFF'
          u = Mo.Von*((e > Mo.db) - (e < -Mo.db));
        case 'P'
          u = S.Kp.*e;
        case 'PI'
          inew = integ + S.Ki.*e*Ts;
          uu = S.Kp.*e + inew;
          free = abs(uu) <= Mo.Vmax;   % anti-windup (conditional integration)
          integ(free) = inew(free);
          u = uu;
        case 'PID'
          inew = integ + S.Ki.*e*Ts;
          uu = S.Kp.*e + inew + S.Kd.*(e-ep)/Ts;
          free = abs(uu) <= Mo.Vmax;
          integ(free) = inew(free);
          u = uu;
      end
      u = max(min(u, Mo.Vmax), -Mo.Vmax);
      ep = e;
    end
    % DC motor: electrical + mechanical (reflected to motor shaft)
    cur = cur + dt*(u - Mo.Ra*cur - Mo.Ke*wm)/Mo.La;
    cur = max(min(cur, Mo.Imax), -Mo.Imax);
    Jt = Mo.Jm + S.Jeff(:,k)./Nr.^2;
    wm = wm + dt*(Mo.Kt*cur - Mo.Bm*wm - S.taud(:,k)./(Nr*Mo.eta))./Jt;
    qm = qm + dt*wm;
    q(:,k+1) = qm./Nr;  Ilog(:,k+1) = cur;  Vlog(:,k+1) = u;
  end
  out.q = q;  out.I = Ilog;  out.V = Vlog;
end

% -------------------------------------------------------------- drawing
function [L, F, Cp] = arm_lines(q, P)
  c1 = cos(q(1)); s1 = sin(q(1));
  u = [c1; s1; 0];  ez = [0; 0; 1];  tv = [-s1; c1; 0];
  p0 = [0;0;0];  p1 = [0;0;P.b1];  p2 = p1 + P.a1*u;
  p3 = p2 + P.a2*(cos(q(2))*u + sin(q(2))*ez);
  q23 = q(2) + q(3);
  p4 = p3 + P.a3*(cos(q23)*u + sin(q23)*ez);
  phi = q23 + q(4);
  d = cos(phi)*u + sin(phi)*ez;
  tip = p4 + P.a4*d;                    % end effector (a4 = 0, so it coincides with the wrist)
  L = [p0 p1 p2 p3 p4 tip];
  A = tip + 0.025*tv;  B = A + 0.05*d;  % parallel gripper fingers ahead of the wrist
  C = tip - 0.025*tv;  D = C + 0.05*d;
  F = [B A C D];
  Cp = tip + 0.025*d;                   % centre of the gripped component
end

function h = draw_arm(q, P, col)
  [L, F, ~] = arm_lines(q, P);
  h(1) = plot3(L(1,:), L(2,:), L(3,:), '-o', 'LineWidth', 3, 'Color', col, ...
               'MarkerFaceColor', col, 'MarkerSize', 6);
  h(2) = plot3(F(1,:), F(2,:), F(3,:), '-', 'LineWidth', 3, 'Color', col);
end

function rack(cx, cy, cz, col)         % shelf / table pad under a component (+ post if raised)
  w = 0.035;  z = cz - 0.02;
  fill3(cx+[-w w w -w], cy+[-w -w w w], z*ones(1,4), col, 'FaceAlpha', 0.4);
  if z > 0.03, plot3([cx cx], [cy cy], [0 z], '-', 'Color', col, 'LineWidth', 2); end
end

function draw_racks(cpick, n)          % storage shelf (blue): two rows, R1..R4
  for k = 1:n
    rack(cpick(1,k), cpick(2,k), cpick(3,k), 'b');
    text(cpick(1,k), cpick(2,k), cpick(3,k)+0.05, sprintf('R%d', k));
  end
end

function draw_bay(cplace, n)           % payload bay (red): ONE square tray, single level, S1..S4
  x0 = min(cplace(1,1:n)) - 0.05;  x1 = max(cplace(1,1:n)) + 0.05;
  y0 = min(cplace(2,1:n)) - 0.05;  y1 = max(cplace(2,1:n)) + 0.05;
  z = cplace(3,1) - 0.02;
  fill3([x0 x1 x1 x0], [y0 y0 y1 y1], z*ones(1,4), 'r', 'FaceAlpha', 0.25);
  plot3([x0 x1 x1 x0 x0], [y0 y0 y1 y1 y0], (z+0.015)*ones(1,5), 'r-', 'LineWidth', 2);   % rim
  for k = 1:n
    plot3(cplace(1,k), cplace(2,k), z+0.001, 'rs', 'MarkerSize', 10);                     % slot marker
    text(cplace(1,k), cplace(2,k), cplace(3,k)+0.05, sprintf('S%d', k));
  end
end

function w = tool_point(c, phi, off)   % end-effector point 'off' metres behind component centre c
  rho = hypot(c(1), c(2));  u = [c(1); c(2)]/rho;       % radial direction
  w = [c(1:2) - off*cos(phi)*u;  c(3) - off*sin(phi)];  % tool axis d = [cos(phi)*u; sin(phi)]
end

% ------------------------------------------------ animation panel callbacks
function stop_anim(src, ~)
  setappdata(ancestor(src, 'figure'), 'stop', true);
end

function play_anim(src, ~)
  fig = ancestor(src, 'figure');
  if getappdata(fig, 'running'), return; end       % press Stop first
  C = getappdata(fig, 'ctx');
  sel = get(C.popup, 'Value');
  nm = C.names{sel};  r = C.R.(nm);
  setappdata(fig, 'stop', false);  setappdata(fig, 'running', true);
  col = C.col(sel,:);
  set(C.hA, 'Color', col, 'MarkerFaceColor', col);  set(C.hF, 'Color', col);
  set(C.hErr, 'XData', C.t(1:50:end), 'YData', max(r.ee(1:50:end), 1e-3));
  tr = zeros(3, 0);
  for k = 1:C.step:numel(C.t)
    if ~ishandle(fig), return; end
    if getappdata(fig, 'stop'), break; end
    [L1, F1, C1] = arm_lines(r.q(:,k), C.P);
    set(C.hA, 'XData', L1(1,:), 'YData', L1(2,:), 'ZData', L1(3,:));
    set(C.hF, 'XData', F1(1,:), 'YData', F1(2,:), 'ZData', F1(3,:));
    for j = 1:C.n
      if k < C.ic(j), cp = C.rack(:,j);
      elseif k < C.io(j), cp = C1;
      else, cp = C.slot(:,j); end
      set(C.hC(j), 'XData', cp(1), 'YData', cp(2), 'ZData', cp(3));
    end
    tr(:, end+1) = L1(:,end);
    set(C.hTr, 'XData', tr(1,:), 'YData', tr(2,:), 'ZData', tr(3,:));
    set(C.hCur, 'XData', C.t(k)*[1 1]);
    set(C.hT, 'String', sprintf('%s control | t = %.2f s | component %d of %d | EE error %.3f mm', ...
        C.label{sel}, C.t(k), min(C.n, 1 + sum(k >= C.io)), C.n, r.ee(k)));
    drawnow; pause(0.002);
  end
  setappdata(fig, 'running', false);
end

%% ======================================================================
%                              MAIN SCRIPT
%% ======================================================================
try, graphics_toolkit('qt'); catch, end   % qt toolkit needed for the control panel
close all; clc;
P  = params();  Mo = motor_params();
rad2d = 180/pi;

%% ---------------------------------------------------- 1. SPECIFICATIONS
fprintf('\n=== 1. DESIGN SPECIFICATIONS ===\n');
fprintf('Configuration : R-R-R-R, base yaw + shoulder + elbow + wrist pitch\n');
fprintf('Links (mm)    : b1=%.0f  a1=%.0f  a2=%.0f  a3=%.0f  a4=%.0f (end effector directly on the wrist)\n', ...
        P.b1*1e3, P.a1*1e3, P.a2*1e3, P.a3*1e3, P.a4*1e3);
fprintf('Reach         : up to %.0f mm radial from the base axis (all links stretched)\n', ...
        (P.a1+P.a2+P.a3+P.a4)*1e3);
fprintf('Joint limits  : ');
for i = 1:4, fprintf('J%d [%4.0f, %4.0f] deg   ', i, P.qmin(i)*rad2d, P.qmax(i)*rad2d); end
fprintf('\nActuators     : 24 V brushed DC, Kt=%.3f N.m/A, Ra=%.1f ohm, La=%.1f mH, gear ratios [%s]\n', ...
        Mo.Kt, Mo.Ra, Mo.La*1e3, num2str(Mo.Nr'));
fprintf('Sensors       : motor-shaft incremental encoder (1024 ppr x4) per joint, homing switch\n');
fprintf('                per joint, gripper opening sensor + fingertip force sensor, current sensing\n');
fprintf('Payload       : %.1f kg (fixed, applied for the whole cycle = conservative)\n', P.mL);

% Tool axis d = [cos(phi)*radial ; sin(phi)] :  phi = -90 deg -> straight down,
% phi = 0 -> horizontal (pointing away from the base), phi = -45 deg -> diagonal.
% STORAGE SHELF (right): two rows - lower row (R1,R2) and upper row (R3,R4). Components 3 and 4
% are grasped from the side / diagonally, so the wrist must rotate before they are set down.
% PAYLOAD BAY (left): a single square tray, one level, with 4 slots at the corners of an 80 mm square;
% every component is set down straight from above.
zh = 0.250;                    % home height
phi = -pi/2;                   % used only for the top-down work-volume check below
NCOMP = 4;            % components handled (>=4 required; use 1 for a quick test run)
cname = {'Avionics board', 'Sensor package', 'Battery unit', 'Comms module'};

% Component RACK on the arm's RIGHT (y<0): 2x2 grid, 80 mm pitch, centres (x;y) in m
rack_xy = [0.24  0.24  0.32  0.32;
          -0.20 -0.28 -0.20 -0.28];
% PAYLOAD-BAY SLOTS on the arm's LEFT (y>0): same 80 mm square footprint, but all on ONE level
slot_xy = [rack_xy(1,:); -rack_xy(2,:)];
zpick     = [0.04 0.04 0.17 0.17];                % component-centre height on the rack (m)
zplace    = [0.04 0.04 0.04 0.04];                % component-centre height in the bay (single level, m)
phi_pick  = deg2rad([-90 -90   0 -45]);           % grasp approach angle
phi_place = deg2rad([-90 -90 -90 -90]);           % insertion approach angle (always straight down)
comp_pick  = [rack_xy; zpick];
comp_place = [slot_xy; zplace];

fprintf('\nComponent assignment (x, y, z of component centre in mm; approach angle in deg)\n');
fprintf('%-3s %-15s %-24s %-7s %-24s %-7s\n', 'No', 'Component', 'Rack position (pick)', 'grasp', 'Payload slot (place)', 'insert');
for k = 1:NCOMP
  fprintf('%-3d %-15s R%d (%4.0f,%5.0f,%4.0f)      %5.0f   S%d (%4.0f,%5.0f,%4.0f)      %5.0f\n', k, cname{k}, ...
          k, comp_pick(1,k)*1e3,  comp_pick(2,k)*1e3,  comp_pick(3,k)*1e3,  phi_pick(k)*180/pi, ...
          k, comp_place(1,k)*1e3, comp_place(2,k)*1e3, comp_place(3,k)*1e3, phi_place(k)*180/pi);
end

% Waypoint matrix (end-effector point at the wrist; the gripped component's centre is 25 mm ahead
% along the tool axis): column 1 = home, then per
% component: above-pick, pick, above-place, place. "above" = backed off 90 mm along the tool axis,
% so lowering/lifting is a straight line along it.
W = [0.300; 0; zh];  phiW = -pi/2;
off = 0.025;  back = 0.09;
for k = 1:NCOMP
  W = [W, tool_point(comp_pick(:,k),  phi_pick(k),  off+back), tool_point(comp_pick(:,k),  phi_pick(k),  off), ...
          tool_point(comp_place(:,k), phi_place(k), off+back), tool_point(comp_place(:,k), phi_place(k), off)];
  phiW = [phiW, phi_pick(k), phi_pick(k), phi_place(k), phi_place(k)];
end
qchk = ik(W(1,:), W(2,:), W(3,:), phiW, P);
if ~all(within_limits(qchk, P))
  error('Waypoint(s) %s unreachable or outside joint limits', mat2str(find(~within_limits(qchk, P))));
end

% Work volume check
[gx, gy, gz] = ndgrid(0.20:0.025:0.40, -0.30:0.05:0.30, 0.06:0.03:0.25);
qg = ik(gx(:)', gy(:)', gz(:)', phi, P);
okg = within_limits(qg, P);
fprintf('Work volume   : x[200,400] y[-300,300] z[60,250] mm -> %.1f %% of %d test points reachable\n', ...
        100*mean(okg), numel(okg));

%% ------------------------------------------------- 2. KINEMATIC ANALYSIS
fprintf('\n=== 2. KINEMATICS ===\n');
fprintf('DH table (frames | theta | b | a | alpha):\n');
fprintf('  1->2 | q1 | %.0f mm | %.0f mm | 90 deg\n', P.b1*1e3, P.a1*1e3);
fprintf('  2->3 | q2 | 0 | %.0f mm | 0\n', P.a2*1e3);
fprintf('  3->4 | q3 | 0 | %.0f mm | 0\n', P.a3*1e3);
fprintf('  4->5 | q4 | 0 | %.0f mm | 0\n', P.a4*1e3);

% FK (DH) vs closed-form vs IK round trip on random configurations
rand('seed', 1);
nt = 200;  maxp = 0;  maxq = 0;  maxc = 0;
for k = 1:nt
  qt = P.qmin + rand(4,1).*(P.qmax - P.qmin);
  T = fk_dh(qt, P);
  pc = fk_closed(qt, P);
  maxc = max(maxc, norm(T(1:3,4) - pc));
  ph = atan2(T(3,1), T(1,1)*cos(qt(1)) + T(2,1)*sin(qt(1)));
  qi = ik(T(1,4), T(2,4), T(3,4), ph, P);
  maxq = max(maxq, max(abs(qi - qt)));
  maxp = max(maxp, norm(fk_dh(qi, P)(1:3,4) - T(1:3,4)));
end
fprintf('Validation on %d random poses:\n', nt);
fprintf('  |FK(DH) - FK(closed form)|   max = %.2e m\n', maxc);
fprintf('  IK->FK position error        max = %.2e m\n', maxp);
fprintf('  IK joint recovery error      max = %.2e rad\n', maxq);

qpick = ik(W(1,3), W(2,3), W(3,3), phiW(3), P);
Tpick = fk_dh(qpick, P);
fprintf('Pick pose joints [deg]: %s\n', num2str(qpick'*rad2d, '%8.2f'));
fprintf('T1-5 at pick pose (m):\n'); disp(Tpick);

%% ------------------------------------------------ 3. PATH PLANNING
fprintf('=== 3. PATH PLANNING ===\n');
dt = 2e-4;                              % plant / reference step
%  segment rows are [type from to T(s)]; per component k (b = 1+4(k-1)):
%    descend (line) | close gripper (dwell) | lift (line) | transfer RIGHT->LEFT (PTP,
%    base swing + wrist re-orientation) | lower (line) | open gripper (dwell) | lift (line) | return (PTP)
segs = [1 1 2 1.5];                     % home -> above first component
for k = 1:NCOMP
  b = 1 + 4*(k-1);
  ap = b+1; pk = b+2; al = b+3; pl = b+4;
  if k < NCOMP, nxt = b+5; Tret = 2.0; else nxt = 1; Tret = 2.0; end
  segs = [segs;
          2 ap pk 0.8;  0 pk pk 0.4;  2 pk ap 0.8;  1 ap al 2.4;
          2 al pl 0.8;  0 pl pl 0.4;  2 pl al 0.8;  1 al nxt Tret];
end
segs = [segs; 0 1 1 0.5];               % settle at home
Tr = make_traj(W, segs, phiW, P, dt, 0.25);
t = Tr.t;  N = numel(t);
if any(isnan(Tr.q(:))) || ~all(within_limits(Tr.q, P))
  error('Path leaves the workspace / joint limits');
end
elz = P.b1 + P.a2*sin(Tr.q(2,:));
fprintf('Wrist (J4) range %.1f to %.1f deg (travel %.0f deg) | lowest elbow height %.0f mm\n', ...
        min(Tr.q(4,:))*rad2d, max(Tr.q(4,:))*rad2d, (max(Tr.q(4,:))-min(Tr.q(4,:)))*rad2d, min(elz)*1e3);
i_close = zeros(1,NCOMP);  i_open = zeros(1,NCOMP);
for k = 1:NCOMP
  s0 = 2 + 8*(k-1);                     % index of the 'descend' segment of component k
  i_close(k) = Tr.iend(s0+1) - round((Tr.iend(s0+1)-Tr.iend(s0))/2);
  i_open(k)  = Tr.iend(s0+5) - round((Tr.iend(s0+5)-Tr.iend(s0+4))/2);
  fprintf('Component %d (%s): grip closes t=%.2f s, released t=%.2f s\n', k, cname{k}, t(i_close(k)), t(i_open(k)));
end
fprintf('Total cycle time for %d components: %.1f s\n', NCOMP, t(end));
fprintf('Peak joint speed [deg/s]: %s\n', num2str(max(abs(Tr.qd),[],2)'*rad2d, '%8.1f'));
fprintf('Peak joint accel [deg/s^2]: %s\n', num2str(max(abs(Tr.qdd),[],2)'*rad2d, '%8.0f'));
pe = fk_closed(Tr.q, P);
vee = sqrt(sum((diff(pe,1,2)/dt).^2, 1));
fprintf('Peak end-effector speed %.2f m/s\n', max(vee));

%% ---------------------------------------- 4. ACTUATOR / DYNAMIC ANALYSIS
fprintf('\n=== 4. DYNAMICS (fixed %.1f kg load) ===\n', P.mL);
sub = 1:10:N;  nc = numel(sub);
tauc = zeros(4,nc);  Jc = zeros(4,nc);  taudc = zeros(4,nc);
for k = 1:nc
  j = sub(k);
  [tk, Mk] = invdyn(Tr.q(:,j), Tr.qd(:,j), Tr.qdd(:,j), P);
  tauc(:,k) = tk;
  Jc(:,k)   = diag(Mk);
  taudc(:,k) = tk - diag(Mk).*Tr.qdd(:,j);   % load seen by each joint (own inertia removed)
end
tc = t(sub);
Jeff = interp1(tc, Jc', t, 'linear')';
taud = interp1(tc, taudc', t, 'linear')';
fprintf('Joint torque [N.m]   peak: %s   rms: %s\n', ...
        num2str(max(abs(tauc),[],2)', '%7.2f'), num2str(sqrt(mean(tauc.^2,2))', '%7.2f'));
Ireq = tauc./(Mo.Nr*Mo.eta*Mo.Kt);
fprintf('Required motor current peak [A]: %s  (limit %.0f A)\n', ...
        num2str(max(abs(Ireq),[],2)', '%7.2f'), Mo.Imax);
wmax = max(abs(Tr.qd),[],2).*Mo.Nr;
fprintf('Peak motor speed [rpm]: %s  (no-load limit %.0f rpm)\n', ...
        num2str(wmax'*30/pi, '%8.0f'), Mo.Vmax/Mo.Ke*30/pi);

%% ----------------------------------------------- 5. CONTROL STRATEGIES
fprintf('\n=== 5. CONTROL ===\n');
Jt   = Mo.Jm + mean(Jeff,2)./Mo.Nr.^2;
wn = 45;  zeta = 0.9;                              % design targets for the PID loop
S.Kp = wn^2*Jt.*Mo.Ra.*Mo.Nr/Mo.Kt;
S.Kd = max((2*zeta*wn*Jt - (Mo.Kt*Mo.Ke/Mo.Ra + Mo.Bm)).*Mo.Ra.*Mo.Nr/Mo.Kt, 0);
S.Ki = 8*S.Kp;
fprintf('Kp [V/rad]: %s\nKi [V/rad/s]: %s\nKd [V.s/rad]: %s\n', ...
        num2str(S.Kp','%9.2f'), num2str(S.Ki','%9.2f'), num2str(S.Kd','%9.3f'));
S.Mo = Mo;  S.dt = dt;  S.t = t;  S.ns = 5;        % controller Ts = 1 ms
S.qref = Tr.q;  S.Jeff = Jeff;  S.taud = taud;

names = {'ONOFF','P','PI','PID'};
R = struct();
tic;
for c = 1:4
  fprintf('Simulating %-6s ...', names{c}); fflush(stdout);
  R.(names{c}) = run_ctrl(names{c}, S);
  fprintf(' done\n');
end
fprintf('(simulation time %.1f s)\n', toc);

% ------ performance metrics
pref = fk_closed(Tr.q, P);
tail = t > t(end) - 0.3;
fprintf('\n%-7s %9s %9s %10s %10s %10s %8s %9s\n', 'Ctrl', 'EE rms', 'EE max', 'err@pick', ...
        'err@place', 'ss err', 'I peak', 'Energy');
fprintf('%-7s %9s %9s %10s %10s %10s %8s %9s\n', '', '[mm]', '[mm]', '[mm]', '[mm]', '[deg]', '[A]', '[J]');
rmsEE = zeros(1,4);
for c = 1:4
  r = R.(names{c});
  pa = fk_closed(r.q, P);
  ee = sqrt(sum((pa - pref).^2, 1))*1e3;
  err = (r.q - Tr.q)*rad2d;
  r.ee = ee;  r.err = err;
  R.(names{c}) = r;
  rmsEE(c) = sqrt(mean(ee.^2));
  Ener = sum(sum(max(r.V.*r.I, 0)))*dt;
  fprintf('%-7s %9.3f %9.3f %10.3f %10.3f %10.4f %8.2f %9.1f\n', names{c}, rmsEE(c), max(ee), ...
          max(ee(i_close)), max(ee(i_open)), mean(mean(abs(err(:,tail)))), max(max(abs(r.I))), Ener);
end
fprintf('(err@pick / err@place = worst case over the %d components)\n', NCOMP);
fprintf('\nEnd-effector error at release, per component [mm]\n%-7s', 'Ctrl');
for k = 1:NCOMP, fprintf(' %14s', cname{k}); end
fprintf('\n');
for c = 1:4
  fprintf('%-7s', names{c});
  for k = 1:NCOMP, fprintf(' %14.3f', R.(names{c}).ee(i_open(k))); end
  fprintf('\n');
end
fprintf('\nPer-joint RMS tracking error [deg]\n');
for c = 1:4
  fprintf('  %-6s %s\n', names{c}, num2str(sqrt(mean(R.(names{c}).err.^2, 2))', '%9.4f'));
end

%% ------------------------------------------------------------- FIGURES
% Fig 1: 3D model, work volume, pick/place poses
kk = min(3, NCOMP);  s0 = 2 + 8*(kk-1);        % snapshots of component 3 (side grasp, wrist rotates)
qs_pick  = Tr.q(:, Tr.iend(s0));
qs_mid   = Tr.q(:, Tr.iend(s0+2) + round((Tr.iend(s0+3)-Tr.iend(s0+2))/2));
qs_place = Tr.q(:, Tr.iend(s0+4));
qs_home  = Tr.q(:, 1);
rand('seed', 2);
qr_ = P.qmin + rand(4,6000).*(P.qmax - P.qmin);
cloud = fk_closed(qr_, P);  cloud = cloud(:, cloud(3,:) > 0);

figure(1); clf; set(gcf, 'Name', '3D model and work volume');
subplot(1,2,1); hold on; grid on;
plot3(cloud(1,:), cloud(2,:), cloud(3,:), '.', 'Color', [0.8 0.8 0.8], 'MarkerSize', 2);
draw_racks(comp_pick, NCOMP);  draw_bay(comp_place, NCOMP);
draw_arm(qs_home,  P, [0.5 0.5 0.5]);
draw_arm(qs_pick,  P, [0 0.3 0.9]);
draw_arm(qs_mid,   P, [0 0.6 0.2]);
draw_arm(qs_place, P, [0.9 0.1 0.1]);
plot3(pref(1,:), pref(2,:), pref(3,:), 'k--', 'LineWidth', 1);
bx = [0.2 0.4]; by = [-0.3 0.3]; bz = [0.06 0.25];
plot3(bx([1 2 2 1 1]), by([1 1 2 2 1]), bz(1)*ones(1,5), 'm-');
plot3(bx([1 2 2 1 1]), by([1 1 2 2 1]), bz(2)*ones(1,5), 'm-');
axis equal; view(40, 25); xlabel('x (m)'); ylabel('y (m)  [-y right, +y left]'); zlabel('z (m)');
title('3D model: pick (blue, right) -> place (red, left)');
subplot(1,2,2); hold on; grid on; axis equal;
plot(cloud(1,:), cloud(2,:), '.', 'Color', [0.8 0.8 0.8], 'MarkerSize', 2);
plot(pref(1,:), pref(2,:), 'k-', 'LineWidth', 1.5);
plot(rack_xy(1,1:NCOMP), rack_xy(2,1:NCOMP), 'bs', 'MarkerSize', 12, 'MarkerFaceColor', 'b');
plot(slot_xy(1,1:NCOMP), slot_xy(2,1:NCOMP), 'rs', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
for k = 1:NCOMP
  text(rack_xy(1,k)+0.012, rack_xy(2,k), sprintf('R%d', k));
  text(slot_xy(1,k)+0.012, slot_xy(2,k), sprintf('S%d %s', k, cname{k}));
end
plot([0.2 0.4 0.4 0.2 0.2], [-0.3 -0.3 0.3 0.3 -0.3], 'm-');
bx0 = min(slot_xy(1,1:NCOMP)) - 0.05;  bx1 = max(slot_xy(1,1:NCOMP)) + 0.05;
by0 = min(slot_xy(2,1:NCOMP)) - 0.05;  by1 = max(slot_xy(2,1:NCOMP)) + 0.05;
plot([bx0 bx1 bx1 bx0 bx0], [by0 by0 by1 by1 by0], 'r-', 'LineWidth', 1.5);
xlabel('x (m)'); ylabel('y (m)'); title('Top view: right (-y) to left (+y)');
legend('reachable', 'EE path', 'pick', 'place', 'base', 'work volume', 'payload bay', 'Location', 'southwest');

% Fig 2: joint motion profiles
figure(2); clf; set(gcf, 'Name', 'Joint motion profiles');
subplot(3,1,1); plot(t, Tr.q*rad2d); grid on; ylabel('q (deg)'); title('Joint motion profiles');
legend('J1','J2','J3','J4','Location','eastoutside');
subplot(3,1,2); plot(t, Tr.qd*rad2d); grid on; ylabel('qdot (deg/s)');
subplot(3,1,3); plot(t, Tr.qdd*rad2d); grid on; ylabel('qddot (deg/s^2)'); xlabel('t (s)');

% Fig 3: Cartesian path
figure(3); clf; set(gcf, 'Name', 'End-effector motion');
subplot(2,1,1); plot(t, pref*1e3); grid on; ylabel('position (mm)');
legend('x','y','z','Location','eastoutside'); title('End-effector motion');
subplot(2,1,2); plot(t(2:end), vee); grid on; ylabel('speed (m/s)'); xlabel('t (s)');

% Fig 4: dynamics
figure(4); clf; set(gcf, 'Name', 'Joint torques and motor current');
subplot(3,1,1); plot(tc, tauc); grid on; ylabel('torque (N.m)');
title('Joint torques (fixed load)'); legend('J1','J2','J3','J4','Location','eastoutside');
subplot(3,1,2); plot(tc, Ireq); hold on; plot(tc([1 end]), [Mo.Imax Mo.Imax], 'k--', tc([1 end]), -[Mo.Imax Mo.Imax], 'k--');
grid on; ylabel('motor current (A)');
subplot(3,1,3); plot(tc, abs(tauc.*Tr.qd(:,sub))); grid on; ylabel('mech. power (W)'); xlabel('t (s)');

% Fig 5: controller comparison
figure(5); clf; set(gcf, 'Name', 'Controller comparison');
sty = {'r-', 'm-', 'g-', 'b-'};
subplot(2,2,1); hold on; grid on;
plot(t, Tr.q(1,:)*rad2d, 'k--', 'LineWidth', 1.5);
for c = 1:4, plot(t, R.(names{c}).q(1,:)*rad2d, sty{c}); end
ylabel('J1 (deg)'); title('Base yaw tracking'); legend('ref', names{:}, 'Location', 'best');
subplot(2,2,2); hold on; grid on;
for c = 1:4, plot(t, R.(names{c}).err(1,:), sty{c}); end
ylabel('J1 error (deg)'); title('Joint 1 error');
subplot(2,2,3); hold on; grid on;
for c = 1:4, semilogy(t, max(R.(names{c}).ee, 1e-4), sty{c}); end
set(gca, 'YScale', 'log'); ylabel('EE error (mm)'); xlabel('t (s)'); title('End-effector position error');
subplot(2,2,4);
bar(rmsEE); set(gca, 'XTickLabel', names, 'YScale', 'log'); grid on;
ylabel('RMS EE error (mm)'); title('RMS error comparison');

% Fig 6: animation panel -- choose a controller in the dropdown, press Run
figure(6); clf;
set(gcf, 'Name', 'Animation panel', 'Position', [100 80 950 720]);
ax = axes('Position', [0.06 0.30 0.90 0.66]); hold on; grid on;
draw_racks(comp_pick, NCOMP);  draw_bay(comp_place, NCOMP);
plot3(pref(1,:), pref(2,:), pref(3,:), 'k:');
axis equal; axis([-0.1 0.6 -0.5 0.5 0 0.6]); view(40, 25);
xlabel('x (m)'); ylabel('y (m)  [-y right, +y left]'); zlabel('z (m)');
[L0, F0, ~] = arm_lines(Tr.q(:,1), P);
hA = plot3(L0(1,:), L0(2,:), L0(3,:), '-o', 'LineWidth', 3, 'Color', [0 0 1], 'MarkerFaceColor', [0 0 1]);
hF = plot3(F0(1,:), F0(2,:), F0(3,:), '-', 'LineWidth', 3, 'Color', [0 0 1]);
hTr = plot3(NaN, NaN, NaN, '-', 'Color', [1 0.4 0], 'LineWidth', 1.2);
ccol = [1 0.6 0; 0 0.7 0.7; 0.6 0.2 0.8; 0.5 0.8 0.2];
hC = zeros(1, NCOMP);
for j = 1:NCOMP
  hC(j) = plot3(comp_pick(1,j), comp_pick(2,j), comp_pick(3,j), 's', 'MarkerSize', 14, ...
                'MarkerFaceColor', ccol(j,:), 'Color', 'k');
end
hT = title('Select a controller and press "Run animation"');

ax2 = axes('Position', [0.08 0.07 0.58 0.17]); hold on; grid on;
hErr = plot(NaN, NaN, 'b-');
hCur = plot([0 0], [1e-3 2e1], 'k-');
set(ax2, 'YScale', 'log'); axis([0 t(end) 1e-3 2e1]);
xlabel('t (s)'); ylabel('EE error (mm)');

C.popup = uicontrol('Style', 'popupmenu', 'String', {'ON-OFF', 'P', 'PI', 'PID'}, 'Value', 4, ...
                    'Units', 'normalized', 'Position', [0.72 0.15 0.22 0.05], 'Callback', @stop_anim);
uicontrol('Style', 'text', 'String', 'Controller:', 'Units', 'normalized', ...
          'Position', [0.72 0.205 0.12 0.035], 'HorizontalAlignment', 'left');
uicontrol('Style', 'pushbutton', 'String', 'Run animation', 'Units', 'normalized', ...
          'Position', [0.72 0.08 0.11 0.05], 'Callback', @play_anim);
uicontrol('Style', 'pushbutton', 'String', 'Stop', 'Units', 'normalized', ...
          'Position', [0.84 0.08 0.10 0.05], 'Callback', @stop_anim);

for c = 1:4                                   % keep only what the animation needs
  Rm.(names{c}).q = R.(names{c}).q;
  Rm.(names{c}).ee = R.(names{c}).ee;
end
C.R = Rm;  C.names = names;  C.label = {'ON-OFF', 'P', 'PI', 'PID'};
C.col = [1 0 0; 0.7 0 0.7; 0 0.6 0; 0 0 1];
C.t = t;  C.P = P;  C.n = NCOMP;  C.ic = i_close;  C.io = i_open;
C.rack = comp_pick;  C.slot = comp_place;  C.step = 250;
C.hA = hA;  C.hF = hF;  C.hC = hC;  C.hTr = hTr;  C.hErr = hErr;  C.hCur = hCur;  C.hT = hT;
setappdata(gcf, 'ctx', C);  setappdata(gcf, 'running', false);  setappdata(gcf, 'stop', false);
fprintf('\nFigure 6: pick a controller in the dropdown, then press "Run animation".\n');
fprintf('Done. Close Figure 6 to finish.\n');
if ~isguirunning()
  waitfor(gcf);        % keeps all windows (and the controller panel) open when run as "octave file.m"
end