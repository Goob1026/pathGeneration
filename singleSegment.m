% Path Generation Baxter
addpath('/home/cats/Documents/oakesk/09_30_2020/pathGen/neededFiles')
clear variables; close all; clc
load('baxWithSledInitPose.mat')
[robot_const, robot_structure] = defineBaxter();
qlimit = [robot_const(1).limit.lower_joint_limit,...
    robot_const(1).limit.upper_joint_limit];
% Specify Lambda (path variable)
N = 100;
lambda = [0:1/N:1];
epsilon_r = 0.1; epsilon_p = 0.1;
q_prime_min = -10*ones(7,1); q_prime_max = 10*ones(7,1);
options = optimoptions('quadprog','Display','off');

%% First Segment
% Specify initial joint angles
q0 = qinit;
% Compute initial pose
[R0T0, P0T0] = fwdkin(robot_const(1).kin,q0); z0 = P0T0(3);
eul0 = rotm2eul(R0T0);
% Specify desired final pose
P0Td = P0T0;%+[-0.2;0.3;0]; 
R0Td = rot([0;0;1],-10*pi/180)*R0T0;
euld = rotm2eul(R0Td);
ep0 = norm(P0Td - P0T0)^2; er0 = 0.5*norm( (sqrtm(R0T0*R0Td') - eye(3)),'fro')^2;
% Compute Path in Task Space
ER0 = R0T0*R0Td'; temp = vrrotmat2vec(ER0); k_hat = [temp(1);temp(2);temp(3)];
theta0 = temp(4); theta_lambda = zeros(1,length(lambda));
dP0T_dlambda = (P0Td - P0T0); % constant in lambda
der_dlambda = zeros(1,length(lambda));
for k = 1:length(lambda)
    theta_lambda(k) = (1 - lambda(k))*theta0;
    der_dlambda(k) = -theta0;
end
% Solve QP Problem and Generate Joint Space Path
q_prime = zeros(7,length(lambda)); q_lambda = zeros(7,length(lambda)); q_lambda(:,1) = q0;
exitflag = zeros(1,length(lambda)); P0T_lambda = zeros(3,length(lambda));
R0T_lambda = zeros(3,3,length(lambda)); P0T_lambda(:,1) = P0T0; R0T_lambda(:,:,1) = R0T0;
eul_lambda = zeros(3,length(lambda)); eul_lambda(:,1) = rotm2eul(R0T0);
rot_dev = zeros(1,length(lambda)); rot_dev(1) = acos(dot(R0T0(:,3),[0;0;1]));
qprev = q0; Ptemp = P0T0; Rtemp = R0T0;
for k = 1:length(lambda)
    [lb,ub] = qprimelimits_full(qlimit,qprev,N,q_prime_max,q_prime_min);
    J = robotjacobian(robot_const(1).kin, qprev);
    
    
    [Aeqrot,beqrot] = equality_constraints_rotation(Rtemp,R0T0',75,J);
    Aeq = [Aeqrot;[J(6,:),0,0]];
    beq = [beqrot;-75*(Ptemp(3)-z0)];

    
    vt = dP0T_dlambda;
    thetax = atan2(Rtemp(3,3),Rtemp(3,2)); 
    thetay = atan2(Rtemp(3,3),Rtemp(3,1));
    vr = der_dlambda(k)*k_hat + 50*[(thetax - pi/2);-(thetay - pi/2);0];
    H = getqp_H(qprev, J, vr, vt, ...
        epsilon_r, epsilon_p);
    f = getqp_f( qprev, epsilon_r, epsilon_p );
    [q_prime_temp,~,exitflag(k)] = quadprog(H,f,[],[],Aeq,beq,lb,ub,[]...
        ,options);
    q_prime_temp = q_prime_temp(1:7);
    % check exit flag - all elements should be 1
    if exitflag(k) ~= 1
        disp('Generation Error')
        return;
    end
    q_prime(:,k) = q_prime_temp; qprev = qprev + (1/N)*q_prime_temp;
    q_lambda(:,k+1) = qprev;   
    [Rtemp, Ptemp] = fwdkin(robot_const(1).kin,qprev);
    P0T_lambda(:,k+1) = Ptemp;  R0T_lambda(:,:,k+1) = Rtemp;
    eul_lambda(:,k+1) = rotm2eul(Rtemp);
    rot_dev(k+1) = acos(dot(Rtemp(:,3),[0;0;1]));
end
% Chop off excess
q_lambda(:,end) = []; P0T_lambda(:,end) = []; R0T_lambda(:,:,end) = []; rot_dev(end) = [];
eul_lambda(:,end) = [];

joint_plots(1,'$\lambda$','q (deg)',lambda,q_lambda*180/pi,qlimit*180/pi)

% Check End Position
qf = q_lambda(:,end);
[R0Tf, P0Tf] = fwdkin(robot_const(1).kin,qf);

R0Td
R0Tf
P0Td
P0Tf


save('SingleSegment.mat','q_lambda','lambda');

% Plot Trajectory in task space (look for error in z-dir)
figure(2)
subplot(2,3,1)
plot(lambda,P0T_lambda(1,:),[lambda(1) lambda(end)],[P0T0(1) P0Td(1)],'LineWidth',2)
xlabel('lambda'); ylabel('x-dir')
subplot(2,3,2)
plot(lambda,P0T_lambda(2,:),[lambda(1) lambda(end)],[P0T0(2) P0Td(2)],'LineWidth',2)
xlabel('lambda'); ylabel('y-dir')
subplot(2,3,3)
plot(lambda,P0T_lambda(3,:),[lambda(1) lambda(end)],[P0T0(3) P0Td(3)],'LineWidth',2)
xlabel('lambda'); ylabel('z-dir')
subplot(2,3,4)
plot(lambda,eul_lambda(3,:)*180/pi,[lambda(1) lambda(end)],[eul0(3) euld(3)]*180/pi,'LineWidth',2)
xlabel('lambda'); ylabel('Roll (deg)')
subplot(2,3,5)
plot(lambda,eul_lambda(2,:)*180/pi,[lambda(1) lambda(end)],[eul0(2) euld(2)]*180/pi,'LineWidth',2)
xlabel('lambda'); ylabel('Pitch (deg)')
subplot(2,3,6)
plot(lambda,eul_lambda(1,:)*180/pi,[lambda(1) lambda(end)],[eul0(1) euld(1)]*180/pi,'LineWidth',2)
xlabel('lambda'); ylabel('Yaw (deg)')

% Plot deviation from contrained orientation
%figure(3)
%plot(lambda,rot_dev*180/pi);title('constrained rotation deviation');ylabel('deg')

rmpath('/home/cats/Documents/oakesk/09_30_2020/pathGen/neededFiles')


