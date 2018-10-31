% set this to be the same as xlims at the top of DDC_mpi.f90
xlims = [0 ; 150];

locfile = fopen('../locs.txt', 'r');
massfile = fopen('../mass.txt', 'r');

shape_locs = fscanf(locfile,'%f', 2);
shape_mass = fscanf(massfile,'%f', 2);

if shape_locs ~= shape_mass
    fprintf('ERROR: locs and mass arrays are different shapes') 
end

locs = fscanf(locfile, '%f');
mass = fscanf(massfile, '%f');

fclose('all');

locs = reshape(locs, shape_locs(1), shape_locs(2));
mass = reshape(mass, shape_mass(1), shape_mass(2));

for i = 1  : shape_locs(2)
    
    figure(1)
    clf
    scatter(locs(:, i), mass(:, i))
    
    if any(locs(:, i) < xlims(1)) || any(locs(:, i) > xlims(2))
        fprintf('ERROR: particle outside boundary in timestep %i \n', i)
    end

    set(gca,'xlim', [xlims(1), xlims(2)], 'ylim', [0, 1]);
    pause(0.02)
end

