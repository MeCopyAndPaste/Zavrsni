function [posD, posA, vel] = beeSimulation_mod3(posD, posA, vel, T, draw, blow)
    
    global casu_pos
    blow_radius_tg  = 4;
    blow_radius_rad = 5;
    % speed
    dist = 0.1;
    % rectangular arena approximation
    arena_len = 15;
    
    % random angular change of bees 
    random_angle = min(max(rand(size(posD))*pi-pi/2,-pi/2),pi/2);
    
    % if left warmer - more towards left; if right warmer - more towards 
    % right.
    % Left always smaller node index.
    left = 1;
    right = 2;
    % if warmer on the left - bees want to go in direction pi - position
    % if warmer on the right - want to go in dir 0 - position
    
    x = posA.*cos(posD);
    y = posA.*sin(posD);
    
   
    tg_DL = 0;
    rad_DL = 0;
    tg_DD = 0;
    rad_DD = 0;
    vel2 = mod(vel, 2*pi);
    %vel 2 je vel u formatu [0 , 2*pi>
    
    %disp('**********************************************************************************************')
    %fprintf('vel_vec     = [%f %f %f %f %f %f %f %f %f %f];\n', vel(1),vel(2),vel(3),vel(4),vel(5),vel(6),vel(7),vel(8),vel(9),vel(10));
    %fprintf('vel2_vec    = [%f %f %f %f %f %f %f %f %f %f];\n', vel2(1),vel2(2),vel2(3),vel2(4),vel2(5),vel2(6),vel2(7),vel2(8),vel2(9),vel2(10));
    %fprintf('x0_vec      = [%f %f %f %f %f %f %f %f %f %f];\n', x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10));
    %fprintf('y0_vec      = [%f %f %f %f %f %f %f %f %f %f];\n', y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8),y(9),y(10));
    
    if blow(1)
        vec_a_DL = atan2( y, x - (-1) * casu_pos );
        vec_a_DL = mod(vec_a_DL, 2*pi);
        vec_a_AL = sqrt(y.^2+ (x - (-1) * casu_pos).^2);
        
        %fprintf('vec_a_DL_vec = [%f %f %f %f %f %f %f %f %f %f];\n', vec_a_DL(1),vec_a_DL(2),vec_a_DL(3),vec_a_DL(4),vec_a_DL(5),vec_a_DL(6),vec_a_DL(7),vec_a_DL(8),vec_a_DL(9),vec_a_DL(10));
        
        tg_DL = vel2- vec_a_DL;
        tg_DL = tg_DL - 2*pi*(tg_DL > pi).*(tg_DL < 2*pi) + 2*pi*(tg_DL < -pi).*(tg_DL>-2*pi);
        tg_DL = vec_a_DL + pi/2 * sign(tg_DL);
        tg_DL = mod(tg_DL, 2*pi);
        
        %fprintf('tg_1L_vec    = [%f %f %f %f %f %f %f %f %f %f];\n', tg_DL(1),tg_DL(2),tg_DL(3),tg_DL(4),tg_DL(5),tg_DL(6),tg_DL(7),tg_DL(8),tg_DL(9),tg_DL(10));
        
        tg_coef_1L  = 0.5*( -(vec_a_AL < blow_radius_tg).*(1-(vec_a_AL./blow_radius_tg))); 
        tg_DL       = vel2 - tg_DL;
        tg_DL       = tg_coef_1L.*(tg_DL - 2*pi*(tg_DL > pi).*(tg_DL < 2*pi) + 2*pi*(tg_DL < -pi).*(tg_DL>-2*pi));
        
        rad_coef_1L =  (vec_a_AL < blow_radius_rad).*(1-(vec_a_AL./blow_radius_rad));
        rad_DL      = rad_coef_1L.*((vec_a_DL-vel2).*(abs(vec_a_DL-vel2) <= pi) + (abs(vec_a_DL-vel2) > pi).*(-1.*sign(vec_a_DL-vel2).*abs(vec_a_DL-vel2)));
        
       
        %fprintf('tg_DL_vec    = [%f %f %f %f %f %f %f %f %f %f];\n', tg_DL(1),tg_DL(2),tg_DL(3),tg_DL(4),tg_DL(5),tg_DL(6),tg_DL(7),tg_DL(8),tg_DL(9),tg_DL(10));
        %fprintf('rad_DL_vec   = [%f %f %f %f %f %f %f %f %f %f];\n', rad_DL(1),rad_DL(2),rad_DL(3),rad_DL(4),rad_DL(5),rad_DL(6),rad_DL(7),rad_DL(8),rad_DL(9),rad_DL(10))
        %fprintf('rad_D2_vec  = [%f %f %f %f %f %f %f %f %f %f];\n', rad_D2(1),rad_D2(2),rad_D2(3),rad_D2(4),rad_D2(5),rad_D2(6),rad_D2(7),rad_D2(8),rad_D2(9),rad_D2(10));
        
    end
    
    if blow(2)
        vec_a_DD = atan2( y, x - casu_pos );
        vec_a_DD = mod(vec_a_DD, 2*pi);
        vec_a_AD = sqrt(y.^2+ (x - casu_pos).^2);
        
        %fprintf('vec_a_DD_vec = [%f %f %f %f %f %f %f %f %f %f];\n', vec_a_DD(1),vec_a_DD(2),vec_a_DD(3),vec_a_DD(4),vec_a_DD(5),vec_a_DD(6),vec_a_DD(7),vec_a_DD(8),vec_a_DD(9),vec_a_DD(10));
        
        tg_DD = vel2- vec_a_DD;
        tg_DD = tg_DD - 2*pi*(tg_DD > pi).*(tg_DD < 2*pi) + 2*pi*(tg_DD< -pi).*(tg_DD>-2*pi);
        tg_DD = vec_a_DD + pi/2 * sign(tg_DD);
        tg_DD = mod(tg_DD, 2*pi);
        
        %fprintf('tg_1D_vec    = [%f %f %f %f %f %f %f %f %f %f];\n', tg_DD(1),tg_DD(2),tg_DD(3),tg_DD(4),tg_DD(5),tg_DD(6),tg_DD(7),tg_DD(8),tg_DD(9),tg_DD(10));
        
        tg_coef_1D  = 0.5*( -(vec_a_AD < blow_radius_tg).*(1-(vec_a_AD./blow_radius_tg))); 
        tg_DD       = vel2 - tg_DD;
        tg_DD       = tg_coef_1D.*(tg_DD - 2*pi*(tg_DD > pi).*(tg_DD < 2*pi) + 2*pi*(tg_DD < -pi).*(tg_DD>-2*pi));
        
        rad_coef_1D =  (vec_a_AD < blow_radius_rad).*(1-(vec_a_AD./blow_radius_rad));
        rad_DD      = rad_coef_1D.*((vec_a_DD-vel2).*(abs(vec_a_DD-vel2) <= pi) + (abs(vec_a_DD-vel2) > pi).*(-1.*sign(vec_a_DD-vel2).*abs(vec_a_DD-vel2)));
        
       
        %fprintf('tg_DD_vec    = [%f %f %f %f %f %f %f %f %f %f];\n', tg_DD(1),tg_DD(2),tg_DD(3),tg_DD(4),tg_DD(5),tg_DD(6),tg_DD(7),tg_DD(8),tg_DD(9),tg_DD(10));
        %fprintf('rad_DD_vec   = [%f %f %f %f %f %f %f %f %f %f];\n', rad_DD(1),rad_DD(2),rad_DD(3),rad_DD(4),rad_DD(5),rad_DD(6),rad_DD(7),rad_DD(8),rad_DD(9),rad_DD(10))
        %fprintf('rad_D2_vec  = [%f %f %f %f %f %f %f %f %f %f];\n', rad_D2(1),rad_D2(2),rad_D2(3),rad_D2(4),rad_D2(5),rad_D2(6),rad_D2(7),rad_D2(8),rad_D2(9),rad_D2(10));
        
    end    
       
 
    
    grad = - (T(left) - T(right) > 0) + (T(left) - T(right) < 0) + 0*(T(left)-T(right) == 0);
    
    phi_temp = atan2(-posA.*sin(posD), grad * casu_pos-posA.*cos(posD));
    phi_temp = phi_temp - vel;
    if phi_temp > pi
        phi_temp = phi_temp - 2 * pi;
    end
    if phi_temp < -pi
        phi_temp = phi_temp + 2 * pi;
    end
%     random_angle = random_angle + phi_temp .* abs(rand(size(random_angle)));
    
    scale = exp(abs(T(left)-T(right))/10 - 1);
    vel = vel + random_angle + (T(left) ~= T(right)) * phi_temp .* ...
        (rand(size(random_angle))).^2 * scale + tg_DL + rad_DL + tg_DD + rad_DD;
    
    if vel < -pi
        vel = vel + 2 * pi;
    end
    if vel > pi
        vel = vel - 2 * pi;
    end
    %fprintf('vel3_vec     = [%f %f %f %f %f %f %f %f %f %f];\n', vel(1),vel(2),vel(3),vel(4),vel(5),vel(6),vel(7),vel(8),vel(9),vel(10));
    
    x = min(max(x + dist * cos(vel),-arena_len/2),arena_len/2);
    y = min(max(y + dist * sin(vel),-arena_len/2),arena_len/2);
    posD = atan2(y,x);
    posA = sqrt(x.^2 + y.^2);
    
    if draw
        scatter(x,y);
        hold on
    %     quiver(x,y,cos(vel),sin(vel));
    %     quiver(x,y,cos(random_angle),sin(random_angle),'--');
    %     quiver(x,y,cos(phi_temp),sin(phi_temp),'--');
        scatter(-casu_pos,0,'x');
        scatter(casu_pos,0,'x');
        hold off
        axis([-10,10,-10,10])
        
    %      hold on;
    %     figure(2)
    %      plot(phi)
    end
end