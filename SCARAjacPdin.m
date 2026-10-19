function Jp=SCARAjacPdin(Q,Qp,L)
% time derivative of the extended jacobian of a 2R robot
l1=L(1); l2=L(2);
g1=L(3); g2=L(4);
Jp=zeros(8,2);
Jp(1,1)=-l1*cos(Q(1))*Qp(1)-l2*cos(Q(1)+Q(2))*(Qp(1)+Qp(2));
Jp(1,2)=-l2*cos(Q(1)+Q(2))*(Qp(1)+Qp(2));
Jp(2,1)=-l1*sin(Q(1))*Qp(1)-l2*sin(Q(1)+Q(2))*(Qp(1)+Qp(2));
Jp(2,2)=-l2*sin(Q(1)+Q(2))*(Qp(1)+Qp(2));
Jp(3,1)=-l1*cos(Q(1))*Qp(1)-g2*cos(Q(1)+Q(2))*(Qp(1)+Qp(2));
Jp(3,2)=-g2*cos(Q(1)+Q(2))*(Qp(1)+Qp(2));
Jp(4,1)=-l1*sin(Q(1))*Qp(1)-g2*sin(Q(1)+Q(2))*(Qp(1)+Qp(2));
Jp(4,2)=-g2*sin(Q(1)+Q(2))*(Qp(1)+Qp(2));
Jp(5,1)=0; Jp(5,2)=0;
Jp(6,1)=-Qp(1)*g1*cos(Q(1)); Jp(6,2)=0;
Jp(7,1)=-Qp(1)*g1*sin(Q(1)); Jp(7,2)=0;
Jp(8,1)=0; Jp(8,2)=0;
end