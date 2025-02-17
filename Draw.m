%% Plot results

% define colors
blue=[0 0 255]/255;
red=[220 20 60]/255;
orange=[255 165 0]/255;
green=[0 205 102]/255;

% start generating pictures
switch settings.model
    
    case 'InvertedPendulum'
        
        figure(4);
        subplot(321)
        plot(time,state_sim(:,1));
        title('p');
        subplot(322)
        plot(time,state_sim(:,2)*180/pi);
        title('\theta');
        subplot(323)
        plot(time,state_sim(:,3));
        title('v');
        subplot(324)
        plot(time,state_sim(:,4)*180/pi);
        title('\omega');
        subplot(3,2,[5 6]);
        title('F');
        stairs(time,controls_MPC(:,1));
        
    case 'ChainofMasses_Lin'
        figure(1);
        subplot(311);
        plot(time,controls_MPC(:,1));
        ylabel('$u_x$','Interpreter','latex');
        subplot(312);
        plot(time,controls_MPC(:,2));
        ylabel('$u_y$','Interpreter','latex');
        subplot(313);
        plot(time,controls_MPC(:,3));
        ylabel('$u_z$','Interpreter','latex');

        n = data.n;
        figure(2);
        plot3([0,state_sim(1,1:n)], [0,state_sim(1,n+1:2*n)], [0,state_sim(1,2*n+1:3*n)],'Color',red,'LineStyle','--');
        hold on;
        grid on;
        plot3([0,state_sim(end,1:n)], [0,state_sim(end,n+1:2*n)],[0,state_sim(end,2*n+1:3*n)],'Color',blue,'LineStyle','-');
        scatter3([0,state_sim(1,1:n)], [0,state_sim(1,n+1:2*n)], [0,state_sim(1,2*n+1:3*n)],10,'MarkerFaceColor','none');
        scatter3([0,state_sim(end,1:n)], [0,state_sim(end,n+1:2*n)],[0,state_sim(end,2*n+1:3*n)],10,'MarkerFaceColor',red);
        xlabel('X[m]');
        ylabel('Y[m]');
        zlabel('Z[m]');
        
        xlim([-0.5 8]);
        ylim([-1 1]);
        zlim([-6 0.5]);
        
   case 'ChainofMasses_NLin'
        figure(1);
        subplot(311);
        plot(time,controls_MPC(:,1));
        ylabel('$u_x$','Interpreter','latex');
        subplot(312);
        plot(time,controls_MPC(:,2));
        ylabel('$u_y$','Interpreter','latex');
        subplot(313);
        plot(time,controls_MPC(:,3));
        ylabel('$u_z$','Interpreter','latex');

        n = data.n;
        figure(2);
        plot3([0,state_sim(1,1:n)], [0,state_sim(1,n+1:2*n)], [0,state_sim(1,2*n+1:3*n)],'Color',red,'LineStyle','--');
        hold on;
        grid on;
        plot3([0,state_sim(end,1:n)], [0,state_sim(end,n+1:2*n)],[0,state_sim(end,2*n+1:3*n)],'Color',blue,'LineStyle','-');
        scatter3([0,state_sim(1,1:n)], [0,state_sim(1,n+1:2*n)], [0,state_sim(1,2*n+1:3*n)],10,'MarkerFaceColor','none');
        scatter3([0,state_sim(end,1:n)], [0,state_sim(end,n+1:2*n)],[0,state_sim(end,2*n+1:3*n)],10,'MarkerFaceColor',red);
        xlabel('X[m]');
        ylabel('Y[m]');
        zlabel('Z[m]');
        
        xlim([-1.2 1.2]);
        ylim([-1.2 1.2]);
                
    case 'TethUAV'
         
        phi_ref = input.od(1,1);
        phi_ref = repmat(phi_ref, size(time));
        theta_ref = input.od(2,1);
        theta_ref = repmat(theta_ref, size(time));
        axes_ref = [];
        axes_lim = [];
        
        figure();      
        subplot(221)
        hold on;
        grid on;
        plot(time(1:end),rad2deg(state_sim(:,1)),'Color',red);
        plot(time(1:end),rad2deg(phi_ref),'k--');
        title('\phi');
        legend('\phi','ref');
        ax = gca; % current axes
        axes_ref = [axes_ref; ax];
        axes_lim = [axes_lim; ax.YLim];
        
        subplot(222)
        hold on;
        grid on;
        plot(time(1:end),rad2deg(state_sim(:,2)),'Color',red);
        title('\phi_{dot}');
        legend('\phi_{dot}')
        ax = gca; % current axes
        axes_ref = [axes_ref; ax];
        axes_lim = [axes_lim; ax.YLim];
         
        subplot(223)
        hold on;
        grid on;
        plot(time(1:end),rad2deg(state_sim(:,3)),'Color',red);
        plot(time(1:end),rad2deg(theta_ref),'k--');
        title('\theta');
        legend('\theta','ref');
        ax = gca; % current axes
        axes_ref = [axes_ref; ax];
        axes_lim = [axes_lim; ax.YLim];
        
        subplot(224)
        hold on;
        grid on;
        plot(time(1:end),rad2deg(state_sim(:,4)),'Color',red);
        title('\theta_{dot}');
        legend('\theta_{dot}');%,'ref');
        ax = gca; % current axes
        axes_ref = [axes_ref; ax];
        axes_lim = [axes_lim; ax.YLim];
        
        % set axes limits, the same for all the plots:
        maxY = max(axes_lim(:,2));
        minY = min(axes_lim(:,1));
        for i = 1 : length(axes_lim)
            cur_ax = axes_ref(i);
            cur_ax.YLim = [minY maxY];
        end
        
        figure()
        subplot(211)
        hold on;
        grid on;
        plot(time(1:end),state_sim(:,5),'Color',red);
        title('f1');
        
        subplot(212)
        hold on;
        grid on;
        plot(time(1:end),state_sim(:,6),'Color',red);
        title('f2');
        
        figure();
        subplot(211)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,1),'Color',red);
        title('df1');
        
        subplot(212)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,2),'Color',red);
        title('df2');
        
        figure();
        grid on;
        plot(time(1:end-1),constraints(:,1),'Color',red);
        title('fL');

        % plot time statistics
        
        figure();
        hold on;
        grid on;
        plot(time(2:end-1)', CPT(2:end, 1)); % cpt, tshooting, tcond, tqp
        title('Time statistics');
        xlabel('[s]')
        ylabel('[ms]')
        
        % plot nonlinear cost fcn terms       
        figure();
        hold on;
        grid on;
        plot(time', rad2deg(state_sim(:,1)+state_sim(:,3)));
        plot(time', rad2deg(pi/2)*ones(size(time)));
        title('Cost fcn: Avoid singularity');
        legend('\phi + \theta', '\pi/2');
        xlabel('[s]')
        ylabel('[deg]')
        
        figure();
        hold on;
        grid on;
        plot(time', rad2deg(state_sim(:,1)), 'k');
        plot(time, rad2deg(phi_ref), 'k--');
        plot(time', rad2deg(state_sim(:,3)), 'r');
        plot(time, rad2deg(theta_ref), 'r--');
        title('Cost fcn: Attitude behavior close to the ground');
        legend('\phi', '\phi_{ref}', '\theta', '\theta_{ref}');
        xlabel('[s]')
        ylabel('[deg]')
        
        figure();
        hold on;
        grid on;
        plot(time', rad2deg(state_sim(:,2)), 'r');
        plot(time', rad2deg(state_sim(:,4)), 'b');
        plot(time', rad2deg(state_sim(:,1)), 'k');
        plot(time, rad2deg(phi_ref), 'k--');
        title('Cost fcn: \phi and \theta velocities close to the ground');
        leg = legend('$\dot{\phi}$', '$\dot{\theta}$', '$\phi$', '$\phi_{ref}$');
        set(leg,'Interpreter','latex');
        xlabel('[s]')
        ylabel('[deg]')
        
    case 'DiM'	

         samples=size(y_sim,1);	

         figure;	
        title('MCA Tracking of perceived signals');	

         subplot(3,2,1);	
        plot(y_sim(:,1),'r');	
        hold on;	
        plot(data.REF(1:samples,1),'k');	
        title('Longitudinal: $\hat{a}_x$','Interpreter','latex');	

         subplot(3,2,2);	
        plot(y_sim(:,2),'r');	
        hold on;	
        plot(data.REF(1:samples,2),'k');	
        title('Lateral: $\hat{a}_y$','Interpreter','latex');	

         subplot(3,2,3);	
        plot(y_sim(:,3),'r');	
        hold on;	
        plot(data.REF(1:samples,3),'k');	
        title('Vertical: $\hat{a}_z$','Interpreter','latex');	

         subplot(3,2,4);	
        plot(y_sim(:,4),'r');	
        hold on;	
        plot(data.REF(1:samples,4),'k');	
        title('Roll: $\hat{\omega}_{\psi}$','Interpreter','latex');	

         subplot(3,2,5);	
        plot(y_sim(:,5),'r');	
        hold on;	
        plot(data.REF(1:samples,5),'k');	
        title('Pitch: $\hat{\omega}_{\theta}$','Interpreter','latex');	

         subplot(3,2,6);	
        plot(y_sim(:,6),'r');	
        hold on;	
        plot(data.REF(1:samples,6),'k');	
        title('Yaw: $\hat{\omega}_{\phi}$','Interpreter','latex');	


         figure;	
        subplot(3,2,1)	
        plot(y_sim(:,13));	
        hold on;	
        plot(y_sim(:,7),'r');	
        plot(zeros(Tf*100,1),'k');	
        title('Longitudinal displacement')	
        lgd=legend('tripod: $p_{x,T}$','hex: $p_{x,H}$','ref: $p_{x,T}$');	
        set(lgd,'Interpreter','latex');	

        subplot(3,2,2)	
        plot(y_sim(:,14));	
        hold on;	
        plot(y_sim(:,8),'r');	
        plot(zeros(Tf*100,1),'k');	
        title('Lateral displacement')	
        lgd=legend('tripod: $p_{y,T}$','hex: $p_{y,H}$','ref: $p_{y,T}$');	
        set(lgd,'Interpreter','latex');	

        subplot(3,2,3)	
        plot(y_sim(:,9),'r');	
        hold on;	
        plot(zeros(Tf*100,1),'k');	
        title('Vertical displacement: $p_{z,H}$','Interpreter','latex');	

        subplot(3,2,4)	
        plot(y_sim(:,20));	
        hold on;	
        plot(y_sim(:,17),'r');	
        plot(zeros(Tf*100,1),'k');	
        title('Yaw')	
        lgd=legend('tripod: $\phi_T$','hex: $\phi_H$','ref');	
        set(lgd,'Interpreter','latex');	

         subplot(3,2,5)	
        plot(y_sim(:,18),'r');	
        hold on;	
        plot(zeros(Tf*100,1),'k');	
        title('Pitch');	
        lgd=legend('hexpod: $\theta_H$','ref');	
        set(lgd,'Interpreter','latex');	

        subplot(3,2,6)	
        plot(y_sim(:,19),'r');	
        hold on;	
        plot(zeros(Tf*100,1),'k');	
        title('Roll');	
        lgd=legend('hexpod: $\psi_H$','ref');	
        set(lgd,'Interpreter','latex');	

        figure;	
        title('Hex actuator constraints')	
        plot(constraints(:,1:6));	
        hold on;	
        plot(1.045*ones(samples,1),':');	
        plot(1.375*ones(samples,1),':');	
        axis([0 mem.iter 1.0 1.4]);	
        title('Hexpod actuator constraints');
        
    case 'TurboEngine'
        
        figure()
        ax1 = subplot(4,1,1);
        hold on
        stairs(time,state_sim(:,3),'Color',blue);
        stairs(time,state_sim(:,4),'--','Color',red);
        ylabel('actuation / %');
        legend('u1', 'u2');
        grid on;

        ax2 = subplot(4,1,2);
        hold on
        plot(time(1:end-1),constraints(:,1),'Color',blue);
        plot(time(1:end-1),data.REF(1)*ones(1,length(time)-1),'k--');
        ylabel('charging pressure / bar');
        legend('output', 'reference')
        grid on;

        ax3 = subplot(4,1,3);
        hold on
        plot(time(1:end-1),constraints(:,2)*60,'Color',blue);
        plot(time(1:end-1),constraints(:,3)*60,'Color',red);
        plot(time(1:end-1),90e3*ones(1,length(time)-1),'--','Color',blue);
        plot(time(1:end-1),180e3*ones(1,length(time)-1),'--','Color',red);
        ylabel('turbocharger speeds / min^{-1}');
        grid on;
        legend('lp', 'hp', 'limit');

        ax4 = subplot(4,1,4);
        hold on
        plot(time,state_sim(:,1),'Color',blue);
        plot(time,state_sim(:,2),'--','Color',red);
        ylabel('x / -');
        legend('x1', 'x2');
        grid on;
        
        xlabel('Time[s]');
        
    case 'Object'
        for i = 1:size(state_sim,1)
            pd_op(:,i) = zyx2R(state_sim(i,4:6))*state_sim(i,7:9)';
            wd_op(:,i) = zyx2R(state_sim(i,4:6))*state_sim(i,10:12)';
        end
        % Position plots
        figure
        grid on
        hold on
        plot(time, p_ref(1:end,1) ,'-','linewidth',2)
        hold on
        plot(time,state_sim(1:end,1),'--','linewidth',2)
        axis([0 time(end) -2 2])
        legend('ref x', 'output x','Interpreter','latex')
        xlabel('t [s]','Interpreter','latex');ylabel('x [m]','Interpreter','latex')
        
        figure
        grid on
        hold on
        plot(time, p_ref(1:end,2) ,'-','linewidth',2)
        hold on
        plot(time,state_sim(1:end,2),'--','linewidth',2)
        axis([0 time(end) -2 2])
        legend('ref y', 'output y','Interpreter','latex')
        xlabel('t [s]','Interpreter','latex');ylabel('y [m]','Interpreter','latex')
        
        figure
        grid on
        hold on
        plot(time, p_ref(1:end,3) ,'-','linewidth',2)
        hold on
        plot(time,state_sim(1:end,3),'--','linewidth',2)
        axis([0 time(end) -2 2])
        legend('ref z', 'output z','Interpreter','latex')
        xlabel('t [s]','Interpreter','latex');ylabel('z [m]','Interpreter','latex')
        
        figure
        grid on
        hold on
        plot3(p_ref(1:end,1),p_ref(1:end,2),p_ref(1:end,3),'-','linewidth',2)
        hold on
        plot3(state_sim(1:end,1),state_sim(1:end,2),state_sim(1:end,3),'--','linewidth',2)
        axis([min(state_sim(1:end,1))-0.5 max(state_sim(1:end,1))+0.5 min(state_sim(1:end,2))-0.5 max(state_sim(1:end,2))+0.5 min(state_sim(1:end,3))-0.5 max(state_sim(1:end,3))+0.5])
        legend('ref traj', 'output traj','Interpreter','latex')
        xlabel('x [m]');ylabel('y [m]');zlabel('z [m]')
        
        % orientation plots
        figure
        grid on
        hold on
        subplot(2,2,1)
        plot(time, o_ref(1:end,1),'-','linewidth',2)
        hold on
        plot(time,state_sim(1:end,4),'--','linewidth',2)
        axis([0 time(end) -2 2])
        legend('Ref Euler Z', 'Output Euler Z')
        xlabel('t [s]','Interpreter','latex');ylabel('Euler angle')
        hold on
        subplot(2,2,2)
        plot(time, o_ref(1:end,2),'-','linewidth',2)
        hold on
        plot(time,state_sim(1:end,5),'--','linewidth',2)
        axis([0 time(end) -2 2])
        legend('Ref Euler Y', 'Output Euler Y')
        xlabel('t [s]');ylabel('Euler angle')
        hold on
        subplot(2,2,3)
        plot(time, o_ref(1:end,3),'-','linewidth',2)
        hold on
        plot(time,state_sim(1:end,6),'--','linewidth',2)
        axis([0 time(end) -2 2])
        legend('Ref Euler X', 'Output Euler X')
        xlabel('t [s]');ylabel('Euler angle')

        % Velocity plots (velocity already rotated so now in W frame)
        figure
        grid on
        hold on
        subplot(2,2,1)
        plot(time, pd_ref(1:end,1), '-','linewidth',2)
        hold on
        plot(time, pd_op(1,1:end),'--','linewidth',2)
        axis([0 time(end) -2 2])
        legend('Ref velocity x', 'Output velocity x')
        xlabel('t [s]','Interpreter','latex');ylabel('vel $[m/s]$','Interpreter','latex')
        hold on
        subplot(2,2,2)
        plot(time, pd_ref(1:end,2), '-','linewidth',2)
        hold on
        plot(time, pd_op(2,1:end),'--','linewidth',2)
        axis([0 time(end) -2 2])
        legend('Ref velocity y', 'Output velocity y')
        xlabel('t [s]');ylabel('vel $[m/s]$','Interpreter','latex')
        hold on
        subplot(2,2,3)
        plot(time, pd_ref(1:end,3), '-','linewidth',2)
        hold on
        plot(time, pd_op(3,1:end),'--','linewidth',2)
        axis([0 time(end) -2 2])
        legend('Ref velocity z', 'Output velocity z')
        xlabel('t [s]');ylabel('vel $[m/s]$','Interpreter','latex')
        
        % Angular Velocity plots
        figure
        grid on
        hold on
        subplot(2,2,1)
        plot(time,w_ref(1,1:end), '-','linewidth',2)
        hold on
        plot(time,wd_op(1,1:end),'--','linewidth',2)
        axis([0 time(end) -2 2])
        legend('Ref ang velocity z', 'Output ang velocity z')
        xlabel('t [s]','Interpreter','latex');ylabel('ang vel $[rad/s]$','Interpreter','latex')
        hold on
        subplot(2,2,2)
        plot(time,w_ref(2,1:end), '-','linewidth',2)
        hold on
        plot(time,wd_op(2,1:end),'--','linewidth',2)
        axis([0 time(end) -2 2])
        legend('Ref ang velocity y', 'Output ang velocity y')
        xlabel('t [s]');ylabel('ang vel $[rad/s]$','Interpreter','latex')
        hold on
        subplot(2,2,3)
        plot(time,w_ref(3,1:end), '-','linewidth',2)
        hold on
        plot(time,wd_op(3,1:end),'--','linewidth',2)
        axis([0 time(end) -2 2])
        legend('Ref ang velocity x', 'Output ang velocity x')
        xlabel('t [s]');ylabel('ang vel $[rad/s]$','Interpreter','latex')
        
        % lambda plots
        figure
        grid minor
        hold on
        subplot(4,4,1)
        % hold on
        plot(time,state_sim(1:end,14),'linewidth',2)
        axis([0 time(end) -2 2])
        xlabel('t [s]','Interpreter','latex');ylabel('$\lambda_1$','Interpreter','latex')
        hold on
        subplot(4,4,2)
        plot(time,state_sim(1:end,15),'linewidth',2)
        axis([0 time(end) -2 2])
        xlabel('t [s]','Interpreter','latex');ylabel('$\lambda_2$','Interpreter','latex')
        hold on
        subplot(4,4,3)
        plot(time,state_sim(1:end,16),'linewidth',2)
        axis([0 time(end) -2 2])
        xlabel('t [s]','Interpreter','latex');ylabel('$\lambda_3$','Interpreter','latex')
        hold on
        subplot(4,4,4)
        plot(time,state_sim(1:end,17),'linewidth',2)
        axis([0 time(end) -2 2])
        xlabel('t [s]','Interpreter','latex');ylabel('$\lambda_4$','Interpreter','latex')
        hold on
        subplot(4,4,5)
        plot(time,state_sim(1:end,18),'linewidth',2)
        axis([0 time(end) -2 2])
        xlabel('t [s]','Interpreter','latex');ylabel('$\lambda_5$','Interpreter','latex')
        hold on
        subplot(4,4,6)
        plot(time,state_sim(1:end,19),'linewidth',2)
        axis([0 time(end) -2 2])
        xlabel('t [s]','Interpreter','latex');ylabel('$\lambda_6$','Interpreter','latex')
        hold on
        subplot(4,4,7)
        plot(time,state_sim(1:end,20),'linewidth',2)
        axis([0 time(end) -2 2])
        xlabel('t [s]','Interpreter','latex');ylabel('$\lambda_7$','Interpreter','latex')
        hold on
        subplot(4,4,8)
        plot(time,state_sim(1:end,21),'linewidth',2)
        axis([0 time(end) -2 2])
        xlabel('t [s]','Interpreter','latex');ylabel('$\lambda_8$','Interpreter','latex')
        hold on
        subplot(4,4,9)
        plot(time,state_sim(1:end,22),'linewidth',2)
        axis([0 time(end) -2 2])
        xlabel('t [s]','Interpreter','latex');ylabel('$\lambda_9$','Interpreter','latex')
        hold on
        subplot(4,4,10)
        plot(time,state_sim(1:end,23),'linewidth',2)
        axis([0 time(end) -2 2])
        xlabel('t [s]','Interpreter','latex');ylabel('$\lambda_{10}$','Interpreter','latex')
        hold on
        subplot(4,4,11)
        plot(time,state_sim(1:end,24),'linewidth',2)
        axis([0 time(end) -2 2])
        xlabel('t [s]','Interpreter','latex');ylabel('$\lambda_{11}$','Interpreter','latex')
        hold on
        subplot(4,4,12)
        plot(time,state_sim(1:end,25),'linewidth',2)
        axis([0 time(end) -2 2])
        xlabel('t [s]','Interpreter','latex');ylabel('$\lambda_{12}$','Interpreter','latex')
        hold on
        subplot(4,4,13)
        plot(time,state_sim(1:end,26),'linewidth',2)
        axis([0 time(end) -2 2])
        xlabel('t [s]','Interpreter','latex');ylabel('$\lambda_{13}$','Interpreter','latex')
        hold on
        subplot(4,4,14)
        plot(time,state_sim(1:end,27),'linewidth',2)
        axis([0 time(end) -2 2])
        xlabel('t [s]','Interpreter','latex');ylabel('$\lambda_{14}$','Interpreter','latex')
        hold on
        subplot(4,4,15)
        plot(time,state_sim(1:end,28),'linewidth',2)
        axis([0 time(end) -2 2])
        xlabel('t [s]','Interpreter','latex');ylabel('$\lambda_{15}$','Interpreter','latex')
        hold on
        subplot(4,4,16)
        % plot(t,pd(1:end-1,1),'-','linewidth',2)
        % hold on
        plot(time,state_sim(1:end,29),'linewidth',2)
        axis([0 time(end) -2 2])
        % legend('Ref velocity x', 'Output velocity x')
        xlabel('t [s]','Interpreter','latex');ylabel('$\lambda_{16}$','Interpreter','latex')

end