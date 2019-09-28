function E3 = E3_test_2(U_field,V_field,x_dist,y_dist)
    len_x = length(V_field(1,:));
    len_y = length(V_field(:,1));
    E3 = zeros(len_x*len_y*8,3);
    for ii = 1:len_x*len_y
%             ii
            V = [U_field(len_y-floor((ii-1)/len_y),mod(ii-1,len_x)+1)...
                V_field(len_y-floor((ii-1)/len_y),mod(ii-1,len_x)+1)];
            if norm(V) ~= 0
                uv_cell_array{ii} = V / norm(V);
            else
                uv_cell_array{ii} = V;
            end
            up_cell_array{ii} = [0;1];
            down_cell_array{ii} = [0;-1];
            right_cell_array{ii} = [1;0];
            left_cell_array{ii} = [-1;0];
            up_right_cell_array{ii} = [x_dist;y_dist]./sqrt(x_dist^2+y_dist^2);
            up_left_cell_array{ii} = [-x_dist;y_dist]./sqrt(x_dist^2+y_dist^2);
            down_right_cell_array{ii} = [x_dist;-y_dist]./sqrt(x_dist^2+y_dist^2);
            down_left_cell_array{ii} = [-x_dist;-y_dist]./sqrt(x_dist^2+y_dist^2);
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

% ii = 1 is bottom left of field, ii = len_x*len_y is top right
    jj = 1;
    for ii = 1:len_x*len_y
            if mod(ii,len_x) ~= 0 %Moving Right
                ind = ii+1;
                cst = right_costs{ii};
                E3(jj,:) = [ii ind cst];
                jj = jj+1;
            end
            if mod(ii,len_x) ~= 1 %Moving to the Left
                ind = ii-1;
                cst = left_costs{ii};
                E3(jj,:) = [ii ind cst];
                jj = jj+1;
            end
            if ii < ((len_x*(len_y - 1)) + 1) %Moving Up
                ind = ii + len_x;
                cst = up_costs{ii};
                E3(jj,:) = [ii ind cst];
                jj = jj+1;
            end
            if ii > len_x %Moving Down
                ind = ii-len_x;
                cst = down_costs{ii};
                E3(jj,:) = [ii ind cst];
                jj = jj+1;
            end
            if ii < ((len_x*(len_y-1)) + 1) &&  mod(ii,len_x) ~= 0 %Moving Up and to the Right
                ind = ii + len_x + 1;
                cst = up_right_costs{ii};
                E3(jj,:) = [ii,ind,cst];
                jj = jj+1;
            end
            if ii < ((len_x*(len_y - 1)) + 1) && mod(ii,len_x) ~= 1 %Moving Up and to the Left
                ind = ii+len_x-1;
                cst = up_left_costs{ii};
                E3(jj,:) = [ii,ind,cst];
                jj = jj+1;
            end
            if ii > len_x && ii >= mod(ii,len_x) ~= 1 %Moving Down and to the Left
                ind = ii-len_x-1;
                cst = down_left_costs{ii};
                E3(jj,:) = [ii,ind,cst];
                jj = jj+1;
            end
            if ii > len_x && mod(ii,len_x) ~= 0 %Moving Down and to the Right
                ind = ii-len_x+1;
                cst = down_right_costs{ii};
                E3(jj,:) = [ii,ind,cst];
                jj = jj+1;
            end
    end
    E3 = E3(1:jj-1,:);
    
end