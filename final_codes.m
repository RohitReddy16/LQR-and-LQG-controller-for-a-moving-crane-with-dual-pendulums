%Two Pendulum Cart 
%By Rohit Reddy Pakhala (rpakhala-119399125) & Sai Teja  

%C. Conditions_for_controlability
clear;
clc;
close all;


defining the variables
syms m1 m2 M g l1 l2 F
%m1, m2 are masses and l1, l2 are lengths respectively
A = [0 1 0 0 0 0; 
    0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0; 
    0 0 -((M*g)+(m1*g))/(M*l1) 0 -g*m2/(M*l1) 0; 
    0 0 0 0 0 1;
    0 0 -m1*g/(M*l2) 0 -((M*g)+(m2*g))/(M*l2) 0];
disp(A)
B = [0;1/M;0;1/(l1*M);0;1/(l2*M)];
disp(B)
%cotrolability matrix is given by
Ct = [B A*B (A^2)*B (A^3)*B (A^4)*B (A^5)*B];
Rank = rank(Ct);
Detetminant = det(Ct);
cc=simplify(det(Ct));
length(A)
ct = Rank - length(A);
%As the result is zero the system is not controllable / Uncontrollable
%So getting a condition for controlability
Ct1 = subs(Ct,l1,l2); %substituting l1=l2
det1 = det(Ct1);
rank(Ct1)
disp('the system is uncontrollable when l1=l2 as det is zero')
disp('the system is controllable when l1 != l2, l1 !=0, and l2 !=0 ')


%D. LQR Controller
As = double(subs(A, {g, l1, l2, m1, m2, M}, {9.8, 20, 10, 100, 100, 1000}));
Bs = double(subs(B, {g, l1, l2, m1, m2, M}, {9.8, 20, 10, 100, 100, 1000}));
cont = ctrb(As,Bs);
if (rank(cont)==6)
    disp("Rank of cont is quals to A Hence the system is controllable")
else
    disp("Rank of cont is not equal to A, system is uncontrollable")
end
% Assumin the Q matrix
Q =[5 0 0 0 0 0;
    0 0 0 0 0 0 ;
    0 0 5000 0 0 0;
    0 0 0 0 0 0 ;
    0 0 0 0 5000 0;
    0 0 0 0 0 0];
R=0.01;
C = eye(6);
D = 0;
%initial state is given below
initial_state = [0;0;20;0;10;0];
%ss is the MATLAB function to find the state space representation of the
%system
sys_1 = ss(As, Bs, C, D);
figure
initial(sys_1,initial_state)
%lqr is in-buit LQR controller function
[K_matrix_gain, S_matrix, P] = lqr(As, Bs, Q, R);
K_matrix_gain;
S_matrix;
P;
%eigen values
eig(S_matrix)
sys_2 = ss(As-(Bs*K_matrix_gain),Bs,C,D);
figure
initial(sys_2,initial_state)

%Non-linear
% x, theta_1 and theta_2 values are defined
%"timespan"
    time = 0:0.1:700; 
%ode45 function is used for definining a differential eqn
    [timeperiod,y1] = ode45(@ydot_func,time,initial_state); 
%plotting output of the function as 2D graph
    plot(timeperiod,y1)
    title('LQR Response to Non-Linear System')
    xlabel('Time')
    ylabel('State Outputs')
    grid on

%E. Observability
%Jacobian linearization of the non linear system of dual pendulum suspended on a crane
%From the problem statement, the output vectors are
% x(t), (theta1(t),theta2(t)), (x(t), theta1(t)),and (x(t),theta1(t),theta2(t)) as C1,C2,C3, and C4 respectively
%because Y = CX+DU, we take C matrix values such that it accounts for the
%state variable
C1 =[1 0 0 0 0 0];
%coressponding to theta1 and theta2
C2 =[0 0 1 0 0 0;
     0 0 0 0 1 0];
%corresponding to x and theta2
C3 =[1 0 0 0 0 0;
     0 0 0 0 1 0];
%corresponding to x, theta1, and theta2
C4 =[1 0 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 0 1 0];
%let observability matrix is defined as O_1, O_2, O_3, and O_4
O_1 = [C1' A'*C1' A'*A'*C1' A'*A'*A'*C1' A'*A'*A'*A'*C1' A'*A'*A'*A'*A'*C1'];
O_2 = [C2' A'*C2' A'*A'*C2' A'*A'*A'*C2' A'*A'*A'*A'*C2' A'*A'*A'*A'*A'*C2'];
O_3 = [C3' A'*C3' A'*A'*C3' A'*A'*A'*C3' A'*A'*A'*A'*C3' A'*A'*A'*A'*A'*C3'];
O_4 = [C4' A'*C4' A'*A'*C4' A'*A'*A'*C4' A'*A'*A'*A'*C4' A'*A'*A'*A'*A'*C4'];
rankArr = [rank(O_1),rank(O_2),rank(O_3),rank(O_4)];
for r = 1:4
    disp("Rank of the observability Matrix: of O_"+r+" is :"+rankArr(r))
    if rankArr(r)==6
        disp('System is observable')
    else 
        disp('System is not observable')
    end
end

%F. Luenberger Observer
%luenberger observer for linearized system
desired_poles = [-1;-3;-2;-4;-6;-8];
%state feedback matrix
k = lqr(As,Bs,Q,R);
% output vectors as per pole placements
L_1 = place(As',C1',desired_poles)';
L_3 = place(As',C3',desired_poles)';
L_4 = place(As',C4',desired_poles)';
Bs_C = [Bs;zeros(size(B))];
% output vector 1
A_C1 = [(As-Bs*k) Bs*k; zeros(size(As)) (As-L_1*C1)];
C_C1 = [C1 zeros(size(C1))];
% ouput vector 3
A_C3 = [(As-Bs*k) Bs*k; zeros(size(As)) (As-L_3*C3)];
C_C3 = [C3 zeros(size(C3))];
%output vector 4
A_C4 = [(As-Bs*k) Bs*k; zeros(size(As)) (As-L_4*C4)];
C_C4 = [C4 zeros(size(C4))];
% intial conditions
initial_state = [1;0;30;2;60;2;zeros(6,1)];
sys1 = ss(A_C1, Bs_C, C_C1,D);
sys3 = ss(A_C3, Bs_C, C_C3,D);
sys4 = ss(A_C4, Bs_C, C_C4,D);
disp("Observer error with x(t) output vector")
figure
initial(sys1,initial_state)
title("Response to initial conditions");
figure
step(sys1)
title("Plotting Step Response for ")

disp("Observer error with x(t), theta2(t)")
figure
initial(sys3,initial_state)
title("Response to initial conditions")
figure
step(sys3)
title("Plotting Step Response")
disp("Observer error with x(t), theta1(t), theta2(t)")
figure
initial(sys4,initial_state)
title("plotting the initial conditions response with Luenberger Observer")
figure
step(sys4)
title("Plotting Step Response")
grid on

%G. LQG Controller
QXU = eye(7);
% process noise
%wgn is in-buit MATLAB to generate gaussian process noise
w = wgn(6,1,4);
% generating measurement noise
v = wgn(1,1,4);
QWV = [w;v]*[w' v'];
%state space representation of the closed loop system with LQG controller
sys5 = ss(As,Bs,C1,D);
sys_LQG = lqg(sys5,QXU,QWV);
initial_state_2 = [1;2;30;3;30;3];
figure
disp("Response at initial conditions")
figure
initial(sys_LQG,initial_state_2)
title("LQG controller for Linear Model")

disp("Response at step LQR control")
figure
step(sys_LQG)
title("LQG controller for Linear Model st step response")
grid on
time_samples = 0:0.1:700;
initial_state_Lqg = [1;2;30;3;30;3;zeros(6,1)];
[t2,y] = ode45(@lqg_func,time_samples,initial_state_Lqg);
%plotting output of the function as 2D graph
plot(time_samples,y)
title('LQR Response to Non-Linear System')
xlabel('Time')
ylabel('State Outputs')
grid on