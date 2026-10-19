function J=SCARAjacdin(Q,L)
% jacobian matrix of the extended 2R
l1=L(1); l2=L(2);
g1=L(3); g2=L(4);
J=zeros(8,2);
J(1,1)=-l1*sin(Q(1))-l2*sin(Q(1)+Q(2));
J(1,2)=-l2*sin(Q(1)+Q(2));
J(2,1)=l1*cos(Q(1))+l2*cos(Q(1)+Q(2));
J(2,2)=l2*cos(Q(1)+Q(2));
J(3,1)=-l1*sin(Q(1))-g2*sin(Q(1)+Q(2));
J(3,2)=-g2*sin(Q(1)+Q(2));
J(4,1)=l1*cos(Q(1))+g2*cos(Q(1)+Q(2));
J(4,2)=g2*cos(Q(1)+Q(2));
J(5,1)=1; J(5,2)=1;
J(6,1)=-g1*sin(Q(1)); J(6,2)=0;
J(7,1)=g1*cos(Q(1)); J(7,2)=0;
J(8,1)=1; J(8,2)=0;
end