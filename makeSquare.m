% Path Generation Baxter

clear variables; close all; clc
addpath([pwd,'/neededFiles'])
load('PoseData.mat')
[robot_const, robot_structure] = defineBaxter();
qlimit = [robot_const(1).limit.lower_joint_limit,...
    robot_const(1).limit.upper_joint_limit];
animate_on = 1;
% Specify Lambda (path variable)
N = 100;
lambda = [0:1/N:1];
epsilon_r = 0.1; epsilon_p = 0.1;
q_prime_min = [-10;-10;-10;-10;-10;-10;-10]; 
q_prime_max = [10;10;10;10;10;10;10];

options = optimoptions('quadprog','Display','off');
qall = zeros(7,404);
lambdaAll = [lambda,1+lambda+1/N,2+2/N+lambda,3+3/N+lambda];

%% First Segment
% Specify initial joint angles
%q0 = qinit(5,:)';
%q0 = [-50;-25;57;47;-104;-7;54]*pi/180;
%q0 = [-24;0;-173;-2;88;-21;85]*pi/180;
q0 = [-51;10;102;50;-29;-20;-67]*pi/180;
% Compute initial pose
[R0T0, P0T0] = fwdkin(robot_const(1).kin,q0); z0 = P0T0(3);
p0 = P0T0;
%p1 = [0.8;0.8;p0(3)];
p1 = [p0(1)-0.3;p0(2)+0.35;p0(3)];
% Specify desired final pose
P0Td = p1; R0Td = R0T0;
ep0 = norm(P0Td - P0T0)^2; er0 = 0.5*norm( (sqrtm(R0T0*R0Td') - eye(3)),'fro')^2;
% Compute Path in Task Space
ER0 = R0T0*R0Td'; temp = vrrotmat2vec(ER0); k_hat = [temp(1);temp(2);temp(3)];
theta0 = temp(4); theta_lambda = zeros(1,length(lambda));
dP0T_dlambda = (P0Td - P0T0); % constant in lambda
der_dlambda = zeros(1,length(lambda));
for k = 1:length(lambda)
    theta_lambda(k) = (1 - lambda(k))*theta0;
    der_dlambda(k) = -2*theta0*sin((theta_lambda(k))/2);
end
% Solve QP Problem and Generate Joint Space Path
q_prime = zeros(7,length(lambda)); q_lambda = zeros(7,length(lambda)); q_lambda(:,1) = q0;
exitflag = zeros(1,length(lambda)); P0T_lambda = zeros(3,length(lambda));
R0T_lambda = zeros(3,3,length(lambda)); P0T_lambda(:,1) = P0T0; R0T_lambda(:,:,1) = R0T0;
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
    rot_dev(k+1) = acos(dot(Rtemp(:,3),[0;0;1]));
end
% Chop off excess
q_lambda(:,end) = []; P0T_lambda(:,end) = []; R0T_lambda(:,:,end) = []; rot_dev(end) = [];

joint_plots(1,'$\lambda$','q (rad)',lambda,q_lambda*180/pi,qlimit*180/pi)

% Check End Position
qf = q_lambda(:,end);
[R0Tf, P0Tf] = fwdkin(robot_const(1).kin,qf);
qall(:,1:101) = q_lambda;

save('Segment1.mat','q_lambda','lambda');

% Plot Trajectory in task space (look for error in z-dir)
figure(2)
subplot(2,2,1)
plot(lambda,P0T_lambda(1,:))
hold on;plot([0 1],[p0(1),p1(1)])
xlabel('lambda')
ylabel('x-dir')
subplot(2,2,2)
plot(lambda,P0T_lambda(2,:))
hold on;plot([0 1],[p0(2),p1(2)])
xlabel('lambda')
ylabel('y-dir')
subplot(2,2,3)
plot(lambda,P0T_lambda(3,:))
hold on;plot([0 1],[p0(3),p1(3)])
xlabel('lambda')
ylabel('z-dir')


% % Plot deviation from contrained orientation
% figure(3)
% plot(lambda,rot_dev)

%% Second Segment
% Specify initial joint angles
q0 = qf;
% Compute initial pose
[R0T0, P0T0] = fwdkin(robot_const(1).kin,q0); z0 = P0T0(3); x0 = P0T0(1); y0 = P0T0(2);
% Specify desired final pose
%p2 = [p1(1);p1(2);0.6];
p2 = [p1(1);p1(2);p1(3)+0.15];
P0Td = p2; R0Td = R0T0;
ep0 = norm(P0Td - P0T0)^2; er0 = 0.5*norm( (sqrtm(R0T0*R0Td') - eye(3)),'fro')^2;
% Compute Path in Task Space
ER0 = R0T0*R0Td'; temp = vrrotmat2vec(ER0); k_hat = [temp(1);temp(2);temp(3)];
theta0 = temp(4); theta_lambda = zeros(1,length(lambda));
dP0T_dlambda = (P0Td - P0T0); % constant in lambda
der_dlambda = zeros(1,length(lambda));
for k = 1:length(lambda)
    theta_lambda(k) = (1 - lambda(k))*theta0;
    der_dlambda(k) = -2*theta0*sin((theta_lambda(k))/2);
end
% Solve QP Problem and Generate Joint Space Path
q_prime = zeros(7,length(lambda)); q_lambda = zeros(7,length(lambda)); q_lambda(:,1) = q0;
exitflag = zeros(1,length(lambda)); P0T_lambda = zeros(3,length(lambda));
R0T_lambda = zeros(3,3,length(lambda)); P0T_lambda(:,1) = P0T0; R0T_lambda(:,:,1) = R0T0;
rot_dev = zeros(1,length(lambda)); rot_dev(1) = acos(dot(R0T0(:,3),[0;0;1]));
qprev = q0; Ptemp = P0T0; Rtemp = R0T0;
for k = 1:length(lambda)
    [lb,ub] = qprimelimits_full(qlimit,qprev,N,q_prime_max,q_prime_min);
    J = robotjacobian(robot_const(1).kin, qprev);
    
    
    [Aeqrot,beqrot] = equality_constraints_rotation(Rtemp,R0T0',75,J);
    Aeq = [Aeqrot;[J(4,:),0,0];[J(5,:),0,0]];
    beq = [beqrot;-75*(Ptemp(1)-x0);-75*(Ptemp(2)-y0)];

    
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
    rot_dev(k+1) = acos(dot(Rtemp(:,3),[0;0;1]));
end
% Chop off excess
q_lambda(:,end) = []; P0T_lambda(:,end) = []; R0T_lambda(:,:,end) = []; rot_dev(end) = [];

joint_plots(11,'$\lambda$','q (rad)',lambda,q_lambda*180/pi,qlimit*180/pi)

% Check End Position
qf = q_lambda(:,end);
[R0Tf, P0Tf] = fwdkin(robot_const(1).kin,qf);
qall(:,102:202) = q_lambda;

save('Segment2.mat','q_lambda','lambda');

% Plot Trajectory in task space (look for error in z-dir)
figure(12)
subplot(2,2,1)
plot(lambda,P0T_lambda(1,:))
hold on;plot([0 1],[p1(1),p2(1)])
xlabel('lambda')
ylabel('x-dir')
subplot(2,2,2)
plot(lambda,P0T_lambda(2,:))
hold on;plot([0 1],[p1(2),p2(2)])
xlabel('lambda')
ylabel('y-dir')
subplot(2,2,3)
plot(lambda,P0T_lambda(3,:))
hold on;plot([0 1],[p1(3),p2(3)])
xlabel('lambda')
ylabel('z-dir')


% Plot deviation from contrained orientation
% figure(13)
% plot(lambda,rot_dev)

%% Third Segment
% Specify initial joint angles
q0 = qf;
% Compute initial pose
[R0T0, P0T0] = fwdkin(robot_const(1).kin,q0); z0 = P0T0(3);
% Specify desired final pose
p3 = [p0(1);p0(2);p2(3)];
P0Td = p3; R0Td = R0T0;
ep0 = norm(P0Td - P0T0)^2; er0 = 0.5*norm( (sqrtm(R0T0*R0Td') - eye(3)),'fro')^2;
% Compute Path in Task Space
ER0 = R0T0*R0Td'; temp = vrrotmat2vec(ER0); k_hat = [temp(1);temp(2);temp(3)];
theta0 = temp(4); theta_lambda = zeros(1,length(lambda));
dP0T_dlambda = (P0Td - P0T0); % constant in lambda
der_dlambda = zeros(1,length(lambda));
for k = 1:length(lambda)
    theta_lambda(k) = (1 - lambda(k))*theta0;
    der_dlambda(k) = -2*theta0*sin((theta_lambda(k))/2);
end
% Solve QP Problem and Generate Joint Space Path
q_prime = zeros(7,length(lambda)); q_lambda = zeros(7,length(lambda)); q_lambda(:,1) = q0;
exitflag = zeros(1,length(lambda)); P0T_lambda = zeros(3,length(lambda));
R0T_lambda = zeros(3,3,length(lambda)); P0T_lambda(:,1) = P0T0; R0T_lambda(:,:,1) = R0T0;
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
    rot_dev(k+1) = acos(dot(Rtemp(:,3),[0;0;1]));
end
% Chop off excess
q_lambda(:,end) = []; P0T_lambda(:,end) = []; R0T_lambda(:,:,end) = []; rot_dev(end) = [];

joint_plots(21,'$\lambda$','q (rad)',lambda,q_lambda*180/pi,qlimit*180/pi)

% Check End Position
qf = q_lambda(:,end);
[R0Tf, P0Tf] = fwdkin(robot_const(1).kin,qf);
qall(:,203:303) = q_lambda;

save('Segment3.mat','q_lambda','lambda');

% Plot Trajectory in task space (look for error in z-dir)
figure(22)
subplot(2,2,1)
plot(lambda,P0T_lambda(1,:))
hold on;plot([0 1],[p2(1),p3(1)])
xlabel('lambda')
ylabel('x-dir')
subplot(2,2,2)
plot(lambda,P0T_lambda(2,:))
hold on;plot([0 1],[p2(2),p3(2)])
xlabel('lambda')
ylabel('y-dir')
subplot(2,2,3)
plot(lambda,P0T_lambda(3,:))
hold on;plot([0 1],[p2(3),p3(3)])
xlabel('lambda')
ylabel('z-dir')


% % Plot deviation from contrained orientation
% figure(23)
% plot(lambda,rot_dev)

%% Fourth Segment
% Specify initial joint angles
q0 = qf;
% Compute initial pose
[R0T0, P0T0] = fwdkin(robot_const(1).kin,q0); z0 = P0T0(3); x0 = P0T0(1); y0 = P0T0(2);
% Specify desired final pose
p4 = p0;
P0Td = p4; R0Td = R0T0;
ep0 = norm(P0Td - P0T0)^2; er0 = 0.5*norm( (sqrtm(R0T0*R0Td') - eye(3)),'fro')^2;
% Compute Path in Task Space
ER0 = R0T0*R0Td'; temp = vrrotmat2vec(ER0); k_hat = [temp(1);temp(2);temp(3)];
theta0 = temp(4); theta_lambda = zeros(1,length(lambda));
dP0T_dlambda = (P0Td - P0T0); % constant in lambda
der_dlambda = zeros(1,length(lambda));
for k = 1:length(lambda)
    theta_lambda(k) = (1 - lambda(k))*theta0;
    der_dlambda(k) = -2*theta0*sin((theta_lambda(k))/2);
end
% Solve QP Problem and Generate Joint Space Path
q_prime = zeros(7,length(lambda)); q_lambda = zeros(7,length(lambda)); q_lambda(:,1) = q0;
exitflag = zeros(1,length(lambda)); P0T_lambda = zeros(3,length(lambda));
R0T_lambda = zeros(3,3,length(lambda)); P0T_lambda(:,1) = P0T0; R0T_lambda(:,:,1) = R0T0;
rot_dev = zeros(1,length(lambda)); rot_dev(1) = acos(dot(R0T0(:,3),[0;0;1]));
qprev = q0; Ptemp = P0T0; Rtemp = R0T0;
for k = 1:length(lambda)
    [lb,ub] = qprimelimits_full(qlimit,qprev,N,q_prime_max,q_prime_min);
    J = robotjacobian(robot_const(1).kin, qprev);
    
    
    [Aeqrot,beqrot] = equality_constraints_rotation(Rtemp,R0T0',75,J);
    Aeq = [Aeqrot;[J(4,:),0,0];[J(5,:),0,0]];
    beq = [beqrot;-75*(Ptemp(1)-x0);-75*(Ptemp(2)-y0)];

    
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
    rot_dev(k+1) = acos(dot(Rtemp(:,3),[0;0;1]));
end
% Chop off excess
q_lambda(:,end) = []; P0T_lambda(:,end) = []; R0T_lambda(:,:,end) = []; rot_dev(end) = [];

joint_plots(31,'$\lambda$','q (rad)',lambda,q_lambda*180/pi,qlimit*180/pi)

% Check End Position
qf = q_lambda(:,end);
[R0Tf, P0Tf] = fwdkin(robot_const(1).kin,qf);
qall(:,304:404) = q_lambda;
save('Segment4.mat','q_lambda','lambda');

% Plot Trajectory in task space (look for error in z-dir)
figure(32)
subplot(2,2,1)
plot(lambda,P0T_lambda(1,:))
hold on;plot([0 1],[p3(1),p4(1)])
xlabel('lambda')
ylabel('x-dir')
subplot(2,2,2)
plot(lambda,P0T_lambda(2,:))
hold on;plot([0 1],[p3(2),p4(2)])
xlabel('lambda')
ylabel('y-dir')
subplot(2,2,3)
plot(lambda,P0T_lambda(3,:))
hold on;plot([0 1],[p3(3),p4(3)])
xlabel('lambda')
ylabel('z-dir')


% Plot deviation from contrained orientation
% figure(33)
% plot(lambda,rot_dev)


joint_plots(41,'$\lambda$','q (rad)',lambdaAll,qall*180/pi,qlimit*180/pi)

load('/home/cats/Documents/oakesk/pathGen/results/Segment1.mat')
for ii = 1:6
subplot(3,3,ii);hold on;plot(lambda,q_lambda(ii,:)*180/pi,'r','LineWidth',2)
end
subplot(3,3,8);hold on;plot(lambda,q_lambda(7,:)*180/pi,'r','LineWidth',2)
load('/home/cats/Documents/oakesk/pathGen/results/Segment2.mat')
for ii = 1:6
subplot(3,3,ii);hold on;plot(lambda+1,q_lambda(ii,:)*180/pi,'r','LineWidth',2)
end
subplot(3,3,8);hold on;plot(lambda+1,q_lambda(7,:)*180/pi,'r','LineWidth',2)
load('/home/cats/Documents/oakesk/pathGen/results/Segment3.mat')
for ii = 1:6
subplot(3,3,ii);hold on;plot(lambda+2,q_lambda(ii,:)*180/pi,'r','LineWidth',2)
end
subplot(3,3,8);hold on;plot(lambda+2,q_lambda(7,:)*180/pi,'r','LineWidth',2)
load('/home/cats/Documents/oakesk/pathGen/results/Segment4.mat')
for ii = 1:6
subplot(3,3,ii);hold on;plot(lambda+3,q_lambda(ii,:)*180/pi,'r','LineWidth',2)
end
subplot(3,3,8);hold on;plot(lambda+3,q_lambda(7,:)*180/pi,'r','LineWidth',2)
legend('new','limit','limit','prev')

rmpath([pwd,'/neededFiles'])


%%
% figure(4);
% baxter = createCombinedRobot(robot_const, robot_structure);
% axis equal;
% axis([-2 2 -2 2 -.85 1]);
% view([109 27]);
% counter = 1;
% if animate_on
%     q = get_angle_structure(baxter);
%     for k = 1:length(lambda)
%         q(1).state = q_lambda(:,k)';
%         q(2).state = [-45 0 0 90 0 0 0]*pi/180;
%         baxter = updateRobot(q, baxter);
%         xlabel('x')
%         ylabel('y')
%         zlabel('z')
%         drawnow
%         M(counter) = getframe(4);
%         counter = counter+1;
%     end
%     
% end

% myVideo = VideoWriter('myfile.avi');
% open(myVideo)
% writeVideo(myVideo, M);
% close(myVideo);