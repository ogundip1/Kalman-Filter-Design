
                       %%%%%%%   MECH 5503 - Adv. Robotics %%%%%%%%%%
                         % Winter 2019, Carleton University %
                     % E-mail: collinsogundipe@cmail.carleton.ca %
                            
%%%% PROBLEM 2 %%%%

clear all
close all

% Actual Position of Boat with time along X-Y coordinates
BoatPose = [(0.1:0.1:10)' (0.2:0.2:20)'];%actual boat position in motion
Transducers = [0.0001  0.0001; 0  22; 4 0; 4 22];%actual transducer positions
sigmatrans_True = 0.005; %% Measurement noise variance
sigmatrans = sigmatrans_True; % Transducers' sonar measurement variance
%sigmatrans = sigmasonar_True*150;%For less optimal model (like 150 times more noise on the model)
sigmadoppler_True = 0.0012;% chosen
sigmadoppler = sigmadoppler_True; %Doppler Velocity variance
%sigmadoppler = sigmadoppler_True/2;%For optimistic model

%simulating measurements
for j = 1:4
    i = 1:100;
    %Range (True/theoretical/desired Range) presumably with no noise
    Range(i,j,:) = sqrt((BoatPose(i,1)-Transducers(j,1)).^2 + (BoatPose(i,2)-Transducers(j,2)).^2);
    
end
%s_measurements
s_measurements = Range; %transducers' measured Range away from the boat by their sonar, w/out noise
s_measurements = s_measurements -sigmatrans_True + 2*sigmatrans_True*rand(100,4);
d_measurements = [linspace(2,2,100)' linspace(0,0,100)'] - sigmadoppler_True + 2*sigmadoppler_True*rand(100,2);%Doppler's velocity boat

%Initialization for least squares
% Linearization of the measurements to be 1-to-1 Linear
x0 = BoatPose';
A = zeros(4,2,100);
l0 = zeros(100,4);    l = zeros(100,4);
w = zeros(100,4);     N = zeros(2,2,100);
U = zeros(2,1,100);   Cx = zeros(2,2,100);
delta = zeros(100,2); xhat = zeros(100,2);
LSQSolutions = zeros(100,2);
P = [(1/sigmatrans^2) 0     0               0
     0     (1/sigmatrans^2) 0    0
     0                0    (1/sigmatrans^2) 0
     0                0    0              (1/sigmatrans^2)];
% Initialization for kalman filter
xk = zeros(100,4); %xk = zeros(4,1);
xkp = zeros(100,4); %xkp = zeros(4,1);
xpredicted = zeros(4,100); %xpredicted = zeros(4,1);
Rk = zeros(4,4,100);%measurement noise %Initialized measurement (GPS) error variance
phi = zeros(4,4);%transition matrix
phi(1,1) = 1; phi(2,2) = 1; phi(3,3) = 1; phi(4,4) = 1;
phi(1,3) = 1;%because the time step is one second
phi(2,4) = 1;%because the time step is one second
sigma_m = 0.91;%% ADJUST this! 0.91 %Standard deviation (Process Noise)
time = linspace(1,100,100)';
param = zeros(100,4); % (Zk - xpredicted)
K_param = zeros(100,4); % K*param => K*(Zk - xpredicted)
% Q => process (Doppler velocity or Acceleration) error variance (noise distribution)
Q = 2*(sigma_m^(2))*[(1/3) 0     (1/2)  0
                   0     (1/3)  0     (1/2)     
                   (1/2)  0     1     0     
                   0      (1/2) 0     1 ];
Pk = zeros(4,4,100);% initial guess for error covariance
Pkpredicted = zeros(4,4,100);
KalmanSolution = zeros(100,4);
predictedSolutions = zeros(100,4);
%looping through all epochs
for i = 1:100
    %filling design matrix and measurement estimates
    %least squares solution for the time-steps
        for j = 1:4
            A(j,1,i) = (BoatPose(i,1)-Transducers(j,1)) ./ ((BoatPose(i,1)-Transducers(j,1)).^2 +(BoatPose(i,2)-Transducers(j,2)).^2).^(1/2);
            A(j,2,i) = (BoatPose(i,2)-Transducers(j,2)) ./ ((BoatPose(i,1)-Transducers(j,1)).^2 +(BoatPose(i,2)-Transducers(j,2)).^2).^(1/2);
            l0(i,j,:) = sqrt((BoatPose(i,1)-Transducers(j,1)).^2 + (BoatPose(i,2)-Transducers(j,2)).^2); 
            l = s_measurements; %noise-induced measured Range
            w(i,j,:) = l0(i,j,:) - l(i,j,:);
        end
        N(:,:,i) = A(:,:,i)'*P*A(:,:,i);  U(:,:,i) = A(:,:,i)'*P*w(i,:,:)';
        delta(i,:) = (-inv(N(:,:,i))*U(:,:,i))';
        xhat(i,1) = x0(1,i) + delta(i,1); xhat(i,2) = x0(2,i) + delta(i,1);
        x0(1,i) = xhat(i,1);  x0(2,i) = xhat(i,2);
        Cx(:,:,i) = inv(N(:,:,i)); Cv(:,:,i) = inv(P) - A(:,:,i)*Cx(:,:,i)*A(:,:,i)';
        LSQSolutions(i,1) = xhat(i,1); LSQSolutions(i,2) = xhat(i,2);

    %kalman filtered solution
    if i == 1
        KalmanSolution(i,1) = xhat(i,1); KalmanSolution(i,2) = xhat(i,2);
        KalmanSolution(i,3) = d_measurements(i,1); KalmanSolution(i,4) = d_measurements(i,2);
        predictedSolutions(i,1) = xhat(i,1); predictedSolutions(i,2) = xhat(i,2);
        predictedSolutions(i,3) = d_measurements(i,1); predictedSolutions(i,4) = d_measurements(i,2);
    end
    if i > 1
        Rk(1,1,i) = Cx(1,1,i); Rk(1,2,i) = Cx(1,2,i); Rk(2,1,i) = Cx(2,1,i); Rk(2,2,i) = Cx(2,2,i);
        Rk(3,3,i) = sigmadoppler^2; Rk(4,4,i) = sigmadoppler^2;
        % Measurement Updates
        Zk(i,1) = xhat(i,1); Zk(i,2) = xhat(i,2);
        Zk(i,3) = d_measurements(i,1); Zk(i,4) = d_measurements(i,2);
        %Prediction
        xkp(i,1) = KalmanSolution(i-1,1); xkp(i,2) = KalmanSolution(i-1,2);
        xkp(i,3) = d_measurements(i-1,1); xkp(i,4) = d_measurements(i-1,2);
        Pkp(1,1,i) = Cx(1,1,i); Pkp(1,2,i) = Cx(1,2,i);
        Pkp(2,1,i) = Cx(2,1,i); Pkp(2,2,i) = Cx(2,2,i);
        Pkp(3,3,i) = sigmadoppler^2; Pkp(4,4,i) = sigmadoppler^2;
        %Prediction
        xpredicted = phi*xkp';
        xpredicted = xpredicted';
        
        predictedSolutions(i,1) = xpredicted(i,1); predictedSolutions(i,2) = xpredicted(i,2);
        predictedSolutions(i,3) = xpredicted(i,3); predictedSolutions(i,4) = xpredicted(i,4);
        
        Pkpredicted = phi.*Pkp.*phi' + Q;
        
        %Kalman Gain
        K(:,:,i) = Pkpredicted(:,:,i)*(Pkpredicted(:,:,i) + Rk(:,:,i))';
        
        %Updating Estimate
        param(i,1) = Zk(i,1) - xpredicted(i,1); param(i,2) = Zk(i,2) - xpredicted(i,2);
        param(i,3) = Zk(i,3) - xpredicted(i,3); param(i,4) = Zk(i,4) - xpredicted(i,4);
        
        K_param(i,:) = (K(:,:,i)*param(i,:)')';
   
        xkp(i,1) = xpredicted(i,1) + K_param(i,1); xkp(i,2) = xpredicted(i,2) + K_param(i,2);
        xkp(i,3) = xpredicted(i,3) + K_param(i,3); xkp(i,4) = xpredicted(i,4) + K_param(i,4);
        %%% KalmanFilter Solutions:
        KalmanSolution(i,1) = xkp(i,1); KalmanSolution(i,2) = xkp(i,2);
        KalmanSolution(i,3) = xkp(i,3); KalmanSolution(i,4) = xkp(i,4);     
    end
end

%RMS Least Squares (LSQ) error computations
LSQRMS_x = sqrt(mean(dot((LSQSolutions(:,1)-BoatPose(:,1)),(LSQSolutions(:,1)-BoatPose(:,1)))));
LSQRMS_y = sqrt(mean(dot((LSQSolutions(:,2)-BoatPose(:,2)),(LSQSolutions(:,2)-BoatPose(:,2)))));
KFRMS_x = sqrt(mean(dot((KalmanSolution(:,1)-BoatPose(:,1)),(KalmanSolution(:,1)-BoatPose(:,1)))));
KFRMS_y = sqrt(mean(dot((KalmanSolution(:,2)-BoatPose(:,2)),(KalmanSolution(:,2)-BoatPose(:,2)))));

%plotting results
figure
subplot(2,1,1)
plot(time, LSQSolutions(:,1))
hold on
plot(time, KalmanSolution(:,1))
plot(time, (predictedSolutions(:,1)))
plot(time, zeros(100,1)+BoatPose(:,1))
legend('LSQ Estimates', 'Kalman Filter', 'Prediction','True Position')
title('X Position Estimates: Measurement Noise Var = 0.2, Process Noise Var = 0.91')
xlabel('time (s)')
ylabel('X position (m)')
ylim([0 12]); 

subplot(2,1,2)
plot(time, LSQSolutions(:,2))
hold on
plot(time, KalmanSolution(:,2))
plot(time, predictedSolutions(:,2))
plot(time, zeros(100,1)+BoatPose(:,2))
legend('LSQ Estimates', 'Kalman Filter', 'Prediction','True Position')
title('Y Position Estimates: Measurement Noise Var = 0.2, Process Noise Var = 0.91')
xlabel('time (s)')
ylabel('Y position (m)')


