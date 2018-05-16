% attributes:
ranges_max = [200, 5, 10, 100, 1, 1, 1, 1, 1, 1, 1, 0.1, 0.1];

for g = 1 : genes
    population(:,g) = rand(population_size, 1) * ...
        (ranges_max(g) - ranges_min(g)) + ranges_min(g);
end

population(:,1) = round(population(:,1));
population(:,4) = round(population(:,4));

ptour = 0.8;
% 1. duration (n * 10 s) 
% 2. intP threshold start initial heating,
% 3. stop initial heating, 
% 4. max duration initial heating
% 5. threshold P heat initial 
% 6. threshold P heat final
% 7. P heat crossover moment
% 8. threshold P cool initial 
% 9. threshold P cool final
% 10. P cool crossover moment
% 11. dT step value


global genes
global ranges_min ranges_max
global ptour pc pmut pcreep

%% initialize

genes = 13;
population_size = 51; % code only supports odd population size
population = zeros(population_size, genes);
new_population = population;
fitness = zeros(population_size, 1);

ranges_min = zeros(1, genes);
ranges_min(1) = 100;
pc = 0.8;
pmut = 1 / genes;
pcreep = 0.8;

iter = 20;

generations = 50;
%% evolution
for gen = 1 : generations 
    disp('--------------------------------------')
    disp(strcat('Generation: ', num2str(gen)))
    %% evaluate
    for ind = 1 : population_size 
        [cor, sub, inc] = Evaluate(population(ind,1:genes), iter);
        fitness(ind) = Objective(cor, sub, inc, iter);
    end
    %% new generation
    % elitism
    [best_fitness, best_ind] = max(fitness);
    
    % debug
    disp(strcat('Fitness of best individual ', num2str(best_fitness)));
    disp('parameters')
    dummy_ind = population(best_ind, :);
    disp(strcat('n = ', num2str(dummy_ind(1))));
    disp(strcat('pIntMin = ', num2str(dummy_ind(2))));
    disp(strcat('pIntMax = ', num2str(dummy_ind(3))));
    disp(strcat('max_ctrl = ', num2str(dummy_ind(4))));
    disp(strcat('start_heat = ', num2str(dummy_ind(5))));
    disp(strcat('stop_heat = ', num2str(dummy_ind(6))));
    disp(strcat('heat_crossover = ', num2str(dummy_ind(7))));
    disp(strcat('start_cool = ', num2str(dummy_ind(8))));
    disp(strcat('stop_cool = ', num2str(dummy_ind(9))));
    disp(strcat('cool_crossover = ', num2str(dummy_ind(10))));
    disp(strcat('rho = ', num2str(dummy_ind(11))));
    disp(strcat('dT_heat = ', num2str(dummy_ind(12))));
    disp(strcat('dT_cool = ', num2str(dummy_ind(13))));
    
    new_population(1,:) = population(best_ind, :);
    % choose other population_size-1 individuals - tournament selection w/
    % crossover and mutation
    for tour = 2 : 2 : (population_size - 1)
        new_ind = Tournament_selection(fitness, population_size);
        individual(1,:) = population(new_ind(1),:);
        individual(2,:) = population(new_ind(2),:);
        % crossover
        individual = Crossover(individual);
        % mutation
        individual = Mutation(individual);
        new_population(tour,:) = individual(1,:);
        new_population(tour+1,:) = individual(2,:);
    end
    population = new_population;
    % number of iterations and max initial heating steps must be integers
    population(:,1) = round(population(:,1));
    population(:,4) = round(population(:,4));
end 

%% functions for code readability

function [cor, sub, inc] = Evaluate(individual, iter)
%     n = individual(1);
%     pIntMin = individual(2);
%     pIntMax = individual(3);
%     max_ctrl = individual(4);
%     start_heat = individual(5);
%     stop_heat = individual(6);
%     heat_crossover = individual(7);
%     start_cool = individual(8);
%     stop_cool = individual(9);
%     cool_crossover = individual(10);
%     rho = individual(11);
%     dT_heat = individual(12);
%     dT_cool = individual(13);

    [cor, sub, inc] = MdsStats(individual, iter);
end

function fitness = Objective(cor, ~, inc, iter)
    % magic number 7 because evaluation on 7 different graph topologies
    fitness = cor/(iter * 7) + 1/(inc + 1);
end

function new_ind = Tournament_selection(fitness, population_size)
    global ptour
    new_ind = [0,0];
    for i = 1 : 2
        ind(1) = round(rand(1) * (population_size) + 0.5);
        ind(2) = round(rand(1) * (population_size) + 0.5);
        p = rand(1);
        [~, winner] = max([fitness(ind(1)), fitness(ind(2))]);
        if p <= ptour
            new_ind(i) = ind(winner);
        else
            % if winner=1, 3-winner=2; else winner=2, 3-winner=1
            new_ind(i) = ind(3 - winner);
        end
    end
end

function individual = Crossover(individual)
    global genes
    global pc
    p = rand(1);
    if p <= pc
        c_point = round(rand(1) * 10 + 0.5);
        ind_temp = [individual(1, 1 : c_point),...
            individual(2, c_point + 1 : genes)];
        individual(2,:) = [individual(2,1:c_point), ...
            individual(1, c_point + 1 : genes)];
        individual(1,:) = ind_temp;
    end
end

function individual = Mutation(individual)
    global genes
    global ranges_min ranges_max
    global pmut pcreep 
    range = ranges_max - ranges_min;
    for i = 1 : 2
        for iGene = 1 : genes
            if rand(1) <= pmut
                if rand(1) <= pcreep
                    % +- 5% around the current value
                    individual(i, iGene) = individual(i, iGene) + ...
                        (rand(1) - 0.5) * range(iGene) * 0.1;
                else
                    % random value in the allowed range
                    individual(i, iGene) = rand(1) * range(iGene) + ...
                        ranges_min(iGene);
                end
            end
        end
    end
end