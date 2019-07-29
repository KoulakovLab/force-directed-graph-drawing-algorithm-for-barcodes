function show_graph2( C, x, y, z, link_col, node_col )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


    %plot3(x,y,z, '.', 'Color', node_col);
    %scatter3(x,y,z, 20, node_col, 'filled');

    [spx, spy, spz] = sphere(10);
    
    cn = sum(C>0,2);
    
    
    
    for i=1:length(x)
        
        r=0.15*(cn(i))^(1/5);
        surf(x(i)+r*spx, y(i)+r*spy, z(i)+r*spz, 'FaceColor', node_col(i,:), 'EdgeColor', 'none', 'FaceLighting', 'phong')
        if i==1
            hold on
        end
        
    end
    
    
    
    
    [I,J] = find(C);
    
    hold on
    for i=1:length(I)
        %if I(i)<=J(i)
        %    continue
        %end
        
        plot3([x(I(i)) x(J(i))], [y(I(i)) y(J(i))], [z(I(i)) z(J(i))], 'Color', link_col)
            
        
    end
    hold off
    axis image

end

