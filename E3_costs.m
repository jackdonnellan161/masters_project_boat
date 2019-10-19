function E3 = E3_test(U_field,V_field,x_dist,y_dist)
    len_x = length(V_field(1,:));
    len_y = length(V_field(:,1));
    E3 = zeros(len_x*len_y*8,3);
    for ii = 1:len_y
        for jj = 1:len_x
            V = [U_field(jj,ii) V_field(jj,ii)];
            if norm(V) ~= 0
                uv_cell_array{jj,ii} = V / norm(V);
            else
                uv_cell_array{jj,ii} = V;
            end
            up_cell_array{jj,ii} = [0;1];
            down_cell_array{jj,ii} = [0;-1];
            right_cell_array{jj,ii} = [1;0];
            left_cell_array{jj,ii} = [-1;0];
            up_right_cell_array{jj,ii} = [1;1]./sqrt(2);
            up_left_cell_array{jj,ii} = [-1;1]./sqrt(2);
            down_right_cell_array{jj,ii} = [1;-1]./sqrt(2);
            down_left_cell_array{jj,ii} = [-1;-1]./sqrt(2);
        end
    end

%     up_costs = cellfun(@(x,y) (-(x*y)+1)*y_dist,uv_cell_array,up_cell_array,...
%         'UniformOutput',false);
%     down_costs = cellfun(@(x,y) (-(x*y)+1)*y_dist,uv_cell_array,down_cell_array,...
%         'UniformOutput',false);
%     left_costs = cellfun(@(x,y) (-(x*y)+1)*x_dist,uv_cell_array,left_cell_array,...
%         'UniformOutput',false);
%     right_costs = cellfun(@(x,y) (-(x*y)+1)*x_dist,uv_cell_array,right_cell_array,...
%         'UniformOutput',false);
%     up_right_costs = cellfun(@(x,y) (-(x*y)+1)*sqrt(x_dist^2+y_dist^2),uv_cell_array,up_right_cell_array,...
%         'UniformOutput',false);
%     up_left_costs = cellfun(@(x,y) (-(x*y)+1)*sqrt(x_dist^2+y_dist^2),uv_cell_array,up_left_cell_array,...
%         'UniformOutput',false);
%     down_right_costs = cellfun(@(x,y) (-(x*y)+1)*sqrt(x_dist^2+y_dist^2),uv_cell_array,down_right_cell_array,...
%         'UniformOutput',false);
%     down_left_costs = cellfun(@(x,y) (-(x*y)+1)*sqrt(x_dist^2+y_dist^2),uv_cell_array,down_left_cell_array,...
%         'UniformOutput',false);

    up_costs = cellfun(@(x,y) acos(x*y)*y_dist,uv_cell_array,up_cell_array,...
        'UniformOutput',false);
    down_costs = cellfun(@(x,y) acos(x*y)*y_dist,uv_cell_array,down_cell_array,...
        'UniformOutput',false);
    left_costs = cellfun(@(x,y) acos(x*y)*x_dist,uv_cell_array,left_cell_array,...
        'UniformOutput',false);
    right_costs = cellfun(@(x,y) acos(x*y)*x_dist,uv_cell_array,right_cell_array,...
        'UniformOutput',false);
    up_right_costs = cellfun(@(x,y) acos(x*y)*sqrt(x_dist^2+y_dist^2),uv_cell_array,up_right_cell_array,...
        'UniformOutput',false);
    up_left_costs = cellfun(@(x,y) acos(x*y)*sqrt(x_dist^2+y_dist^2),uv_cell_array,up_left_cell_array,...
        'UniformOutput',false);
    down_right_costs = cellfun(@(x,y) acos(x*y)*sqrt(x_dist^2+y_dist^2),uv_cell_array,down_right_cell_array,...
        'UniformOutput',false);
    down_left_costs = cellfun(@(x,y) acos(x*y)*sqrt(x_dist^2+y_dist^2),uv_cell_array,down_left_cell_array,...
        'UniformOutput',false);


    jj = 1;
    for ii = 1:len_x*len_y
        [x_cost,y_cost] = ind2sub([len_x len_y],...
            ii);
            if x_cost ~= len_x %Moving to the Right
                ind = sub2ind([len_x len_y],x_cost+1,y_cost);
                cst = right_costs{x_cost,y_cost};
                E3(jj,:) = [ii ind cst];
                jj = jj+1;
            end
            if x_cost ~= 1 %Moving to the Left
                ind = sub2ind([len_x len_y],x_cost-1,y_cost);
                cst = left_costs{x_cost,y_cost};
                E3(jj,:) = [ii ind cst];
                jj = jj+1;
            end

            if y_cost ~= len_y %Moving Up
                ind = sub2ind([len_x len_y],x_cost,y_cost+1);
                cst = up_costs{x_cost,y_cost};
                E3(jj,:) = [ii ind cst];
                jj = jj+1;
            end
            if y_cost ~= 1 %Moving Down
                ind = sub2ind([len_x len_y],x_cost,y_cost-1);
                cst = down_costs{x_cost,y_cost};
                E3(jj,:) = [ii ind cst];
                jj = jj+1;
            end
            if (x_cost ~= len_x && y_cost ~= len_y) %Moving Up and to the Right
                ind = sub2ind([len_x len_y],x_cost+1,y_cost+1);
                cst = up_right_costs{x_cost,y_cost};
                E3(jj,:) = [ii,ind,cst];
                jj = jj+1;
            end
            if (x_cost ~= 1 && y_cost ~= len_y) %Moving Up and to the Left
                ind = sub2ind([len_x len_y],x_cost-1,y_cost+1);
                cst = up_left_costs{x_cost,y_cost};
                E3(jj,:) = [ii,ind,cst];
                jj = jj+1;
            end
            if (x_cost ~= 1 && y_cost ~= 1) %Moving Down and to the Left
                ind = sub2ind([len_x len_y],x_cost-1,y_cost-1);
                cst = down_left_costs{x_cost,y_cost};
                E3(jj,:) = [ii,ind,cst];
                jj = jj+1;
            end
            if (x_cost ~= len_x && y_cost ~= 1) %Moving Down and to the Right
                ind = sub2ind([len_x len_y],x_cost+1,y_cost-1);
                cst = down_right_costs{x_cost,y_cost};
                E3(jj,:) = [ii,ind,cst];
                jj = jj+1;
            end
    end
    E3 = E3(1:jj-1,:);
    
end