function PlotScara(Q,L,col,fig,ll)
q1=Q(1); q2=Q(2);
l1=L(1); l2=L(2);
x1=l1*cos(q1);
y1=l1*sin(q1);
x2=x1+l2*cos(q1+q2);
y2=y1+l2*sin(q1+q2);
figure(fig); grid on;
set(ll,'XData', [0 x1 x2], 'YData', [0 y1 y2], 'color', col);
axis equal;
app=(l1+l2)*1.1;
xlim([-app app]);
ylim([-app app]);
drawnow;
end