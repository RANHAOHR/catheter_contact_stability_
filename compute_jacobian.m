function [J] = compute_jacobian(theta)


wz = [0,0,1]';
wy = [0,-1,0];
wx = [1,0,0]';

q6 = [33,-13,14.7]';
q5 = [24,-13,14.7]';
q4 = [12,-13,14.7]';
q3 = [0,-13,14.7]';
q2 = [0,-13,14.7]';
q1 = [0,-13,14.7]';

%e6 = compute_ew(wx,theta(6));
e1 =compute_exponential_map(wz,theta(1),q1);
e2 =compute_exponential_map(wy,theta(2),q2);
e3 =compute_exponential_map(-wz,theta(3),q3);
e4 =compute_exponential_map(-wz,theta(4),q4);
e5 =compute_exponential_map(-wz,theta(5),q5);

q1_ = [0,-13,0]';
q2_ = [0,-13,14.7]';
q3_ = e1*e2*[0,-13,14.7,1]';

q4_ = e1*e2*e3*[12,-13,14.7,1]';
q5_ = e1*e2*e3*e4*[24,-13,14.7,1]';
q6_ = e1*e2*e3*e4*e5*[33,-13,14.7,1]';

w1_ = [0,0,1]';
w2_ = e1*[0,-1,0,0]';
w3_ = e1*e2*[0,0,-1,0]';
w4_ = e1*e2*e3*[0,0,-1,0]';
w5_ = e1*e2*e3*e4*[0,0,-1,0]';
w6_ = e1*e2*e3*e4*e5*[1,0,0,0]';

q3_ = q3_(1:3);
q4_ = q4_(1:3);
q5_ = q5_(1:3);
q6_ = q6_(1:3);
w2_ = w2_(1:3);
w3_ = w3_(1:3);
w4_ = w4_(1:3);
w5_ = w5_(1:3);
w6_ = w6_(1:3);

j1 = [-cross(w1_,q1_);w1_];
j2 = [-cross(w2_,q2_);w2_];
j3 = [-cross(w3_,q3_);w3_];
j4 = [-cross(w4_,q4_);w4_];
j5 = [-cross(w5_,q5_);w5_];
j6 = [-cross(w6_,q6_);w6_];

J = [j1,j2,j3,j4,j5,j6];



