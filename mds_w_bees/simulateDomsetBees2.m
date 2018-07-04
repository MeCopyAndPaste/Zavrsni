function [vec, pInt] = simulateDomsetBees2(N,varargin)

%% initialize
global casu_pos
numvarargs = length(varargin);

if numvarargs > 3
    error('simulateDomset requires at most 2 optional inputs');
end

optargs = {200 1};
optargs(1:numvarargs) = varargin;
[n, rho] = optargs{:};

k = 50;
T = 28 * ones(size(N));
pBuff = zeros([size(N), k]);
P = zeros(size(N));
pInt = P;

nBees = 10;
posD = randn(sum(sum(N))/2,nBees);
posA = randn(sum(sum(N))/2,nBees);
vel = randn(sum(sum(N))/2,nBees);

ctrl = zeros(length(N),1);

heat = ctrl';
heat_float = ctrl';
cool = ctrl';
cool_float = ctrl';

stop_heat = 0.7;
start_heat = 0.1;
stop_cool = 0.5;
start_cool = 0.2;

step_cool_modifier = ones(size(N));

blowbuff = [];
blow_start = [];
blow = zeros(size(N));

minPbuff = [];

draw = 0;

%% Prepare figure for bee simulation
if draw
    for ifig = 1 : sum(sum(N))/2
        f = figure(ifig);
        grid on
        axis([-10,10,-10,10])
        str = 't = 0';
        annotation('textbox',[0.7,0,1,0.9],'String',str,'FitBoxToText','on','Tag','time_tag');

        mTextBox1 = uicontrol('style','text','Tag','textbox1');
        mString1 = sprintf('T = %.2f', T(1));
        set(mTextBox1,'String',mString1);
        mTextBoxPosition1 = get(mTextBox1,'Position');
        set(mTextBox1,'Position',[150,330,100,15]);

        mTextBox2 = uicontrol('style','text','Tag','textbox2');
        mString2 = sprintf('T = %.2f', T(2));
        set(mTextBox2,'String',mString2);
        mTextBoxPosition2 = get(mTextBox2,'Position');
        set(mTextBox2,'Position',[350,330,100,15]);
    end
end
%% (n * k / 10 seconds)
for i = 1 : n
    iPair = 0;
    
    for iSmaller = 1 : length(N)
        for iLarger = iSmaller + 1 : length(N)
            if N(iSmaller,iLarger) == 1
                iPair = iPair + 1;
                Tsubarena = [T(iSmaller),T(iLarger)];
                blowsubarena = [blow(iSmaller),blow(iLarger)];
                %% move bees k times
                for iTime = 1 : k
                    if draw
                        fi = figure(iPair);
                    end
                    [posD(iPair,:),posA(iPair,:),vel(iPair,:)] = beeSimulation2...
                        (posD(iPair,:),...
                        posA(iPair,:),...
                        vel(iPair,:),...
                        Tsubarena, ...
                        draw, blowsubarena);
                    %% update plot info
                    if draw
                        title(num2str(iPair));
                        delete(findall(fi,'Tag','time_tag'));
                        str = ['t = ',int2str(i)];
                        annotation('textbox',[0.7,0,1,0.9],'String',str,'FitBoxToText','on','Tag','time_tag');
                        mString1 = sprintf('T = %.2f', Tsubarena(1));
                        m = findall(fi,'Tag','textbox1');
                        set(m,'String',mString1);
                        mString2 = sprintf('T = %.2f', Tsubarena(2));
                        m = findall(fi,'Tag','textbox2');
                        set(m,'String',mString2);
                        axis([-10,10,-10,10])
                        pause(0.0005);
                        
                        if ~ishghandle(f)
                            close all
                            break
                        end
                    end
                    pBuff(iSmaller, iLarger, iTime) = ...
                        sum((posA(iPair,:).*cos(posD(iPair,:)) + casu_pos).^2 +...
                        (posA(iPair,:).*sin(posD(iPair,:))).^2 < 4) / nBees;
                    pBuff(iLarger, iSmaller, iTime) = ...
                        sum((posA(iPair,:).*cos(posD(iPair,:)) - casu_pos).^2 +...
                        (posA(iPair,:).*sin(posD(iPair,:))).^2 < 4) / nBees;
                end
%                 P(iSmaller,iLarger) = ...
%                     sum((posA(iPair,:).*cos(posD(iPair,:)) + casu_pos).^2 +...
%                     (posA(iPair,:).*sin(posD(iPair,:))).^2 < 4) / nBees;
%                 P(iLarger,iSmaller) = ...
%                     sum((posA(iPair,:).*cos(posD(iPair,:)) - casu_pos).^2 +...
%                     (posA(iPair,:).*sin(posD(iPair,:))).^2 < 4) / nBees;
                  P(iSmaller, iLarger) = mean(pBuff(iSmaller, iLarger, :), 3);
                  P(iLarger, iSmaller) = mean(pBuff(iLarger, iSmaller, :), 3);
              end
        end
    end
    pInt = pInt + P;
%     P
    ctrl = ctrl + ((sum(pInt,2) > 2) & (ctrl < 10));
    
    arenaCnt = sum(N,2);
    for iNode = 1 : length(N)
        maxP(iNode) = max(P(iNode, N(iNode, :)==1));
        avgP(iNode) = mean(P(iNode, N(iNode, :)==1));
        minP(iNode) = min(P(iNode, N(iNode, :)==1));
    end
    
% lines 140 to 151 for testing purposes only    
%     if i > 80
%        minP = 0.5*ones(length(N));
%        minP = minP(1, :);
%     end
%     if i > 110
%        minP = 0*ones(length(N));
%        minP = minP(1, :);
%     end 
%     if i > 150
%        minP = 0.5*ones(length(N));
%        minP = minP(1, :);
%     end

    progress_smooth_heat = 1 - exp(0.17 * (1 - 1 / (1 - min(1, 3 * i / n ))));
    progress_smooth_cool = 1 - exp(0.85 * (1 - 1 / (1 - min(1, 2 * i / n))));
    %progress_smooth_blow = 1 - exp(0.5 * (1 - 1 / (1 - min(1, 2 * i / n ))));
    
    
    scaling_heat = (1-progress_smooth_heat) * start_heat + stop_heat * progress_smooth_heat;
    scaling_cool = (1-progress_smooth_cool) * start_cool + stop_cool * progress_smooth_cool;
    %scaling_blow = (1-progress_smooth_blow) * start_blow + stop_blow * progress_smooth_blow;
    
    
    
    minPbuff  = [minPbuff; minP];
    blowbuff  = [blowbuff; transpose(blow(:,1))];
    blow_start= [blow_start; 0*minP];
    %    t > 5 min     t < 15 min  after simulation start
    if ( i > 60 )
        for blowi = 1:length(N)
            if (blowbuff( i , blowi) == 1) && (blowbuff( i-1 , blowi) == 0)
                blow_start( i, blowi) = 1;
            else    
                blow_start( i, blowi) = 0;
            end
        end
    end
    if ( i >= 60 ) && ( i <= 180 ) 
        for blowi = 1:length(N)
            minPcomp = max(minPbuff((i-12):i,blowi)); 
            blow(blowi,1)= (minPcomp < 0.2) || ( sum(blow_start((i-18):i,blowi)) ~= 0);
            blow_start(i, blowi) = (minPcomp < 0.2) && ( sum(blow_start((i-18):i,blowi)) ~= 0);
        end
    else
        blow = zeros(size(N));
    end
    
    
    cool_float = cool_float * (1-rho) + (maxP < scaling_cool .* (ctrl > 0)') * rho;
    cool = cool_float > 0.5;
    
    heat_float = heat_float * (1-rho) + ((avgP > scaling_heat) .* ...
        (ctrl > 0)' & (cool == 0)) * rho;
    heat = heat_float > 0.5;
   

    T(:,1) = min(36, T(:,1) + min(0.5, 0.05 * arenaCnt .* (heat')));
    
    step_cool_modifier = step_cool_modifier + 2 * blow;    
    T(:,1) = max(26, T(:,1) - 0.03.*step_cool_modifier(:,1).*(cool'));% - arenaCnt == zeros(size(cool))));
    
    step_cool_modifier = ones(size(N));
    vec(i,:) = T(:,1);
    
    for c = 1 : length(N)
        T(c,:) = T(c,1);
        blow(c,:) = blow(c,1);
    end
    
end
% ploton = 1;
%      if ploton
%      figure(99);
%      clf;
%      hold on;
%          for znj = 1 :length(N)
%                      subplot(length(N),1, znj);
%                      plot(1:n,blowbuff(1:n,znj), 'b',1:i,minPbuff(1:n,znj), 'r');
%  
%          end
%      hold off;   
%      end
% disp("Press any key to continue");
% pause;
pInt = pInt / n;

end
