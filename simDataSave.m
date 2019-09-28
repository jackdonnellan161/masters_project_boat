test_name = 'divide_by_zero_fix_check';
mkdir(test_name)
cd(test_name)
figure
contour(xs,ys,Psi_in,10)
hold on
scatter(global_points(:,1),global_points(:,2),'r')
scatter(waypoints(:,1),waypoints(:,2),'b')
scatter(global_points(1,1),global_points(1,2),'go')
scatter(global_points(end,1),global_points(end,2),'gx')
axis('equal')
scatter(simout_x_boat(:,1),simout_x_boat(:,2),'c*')
legend('Flow Streamlines','Planned Path','Waypoints','Start','End',...
    'Actual Path')
title(strcat('Path Overview: ',test_name))
saveas(gcf,strcat(test_name,'_path overview.png'))

filename = strcat(test_name,'.mat');
save(filename,'simout_abs_T_sum','simout_a_boat','simout_heading',...
    'simout_rot_accels','simout_rot_rates','simout_v_boat','simout_x_boat')