clear variables

% 2 mass transfer modes (i.e., SPaM and PLSS)
nXfer = 2;

DDCTime = '../PDDC_times.txt';
DDCErr  = '../PDDC_Error.txt';

%% read the files in and put them in sensibly-shaped arrays

timefile = fopen(DDCTime, 'r');
errfile  = fopen(DDCErr,  'r');

% get the shape of the arrays
shape_time = fgetl(timefile);
shape_time = str2num(shape_time);
shape_err  = fgetl(errfile);
shape_err  = str2num(shape_err);

if sum(shape_time - shape_err) ~= 0
    fprintf('ERROR: time file is not the same length as error file \n')
    return
else
    shapemat = shape_time;
end

nEns = shapemat(1);
nNp  = shapemat(2);

% these matrices are dimension num_ens x num_particle_ens x nXfer (= 2)
TimeMat = zeros(nEns, nNp, nXfer);
ErrMat  = zeros(nEns, nNp, nXfer);

% get the particle number list
Np_list_time = fgetl(timefile);
Np_list_time = str2num(Np_list_time);
Np_list_err  = fgetl(errfile);
Np_list_err  = str2num(Np_list_err);

if sum(Np_list_time - Np_list_err) ~= 0
    fprintf('ERROR: Np list differs for time and error files \n')
    return
else
    Np_list = Np_list_time;
end

% get the linear vectors of time and error
tmpTime = fscanf(timefile, '%f');
tmpErr  = fscanf(errfile,  '%f');

% reshape the vectors into 3D arrays
ii = 1;
for i = 1 : nNp
    for j = 1 : nXfer
        for k = 1 : nEns
            TimeMat(k, i, j) = tmpTime(ii);
            ErrMat(k, i, j)  = tmpErr(ii);
            ii = ii + 1;
        end
    end
end

fclose('all');

%% average the ensembles

% row 1 is SPaM and row 2 is PLSS
DDCTime = zeros(nXfer, length(Np_list));
DDCErr  = zeros(nXfer, length(Np_list));

for i = 1 : nXfer
    DDCTime(i, :) = mean(TimeMat(:, :, i), 1);
    DDCErr(i, :)  = mean(ErrMat(:, :, i),  1);
end

%% plotting

% Time plot (linear)
figure(1)
clf
plot(Np_list, DDCTime(1, :))
hold on
plot(Np_list, DDCTime(2, :))
title('Run Time (linear)', 'Interpreter', 'latex', 'FontSize', 20)
xlabel('Number of Particles ($N$)', 'Interpreter', 'latex', 'FontSize', 18)
ylabel('Run Time (s)', 'Interpreter', 'latex', 'FontSize', 18)
legend({'SPaM', 'PLSS'}, 'Interpreter', 'latex', 'FontSize', 16,'Location', 'best')

% Time plot (loglog)
figure(2)
clf
loglog(Np_list, DDCTime(1, :))
hold on
loglog(Np_list, DDCTime(2, :))
title('Run Time (loglog)', 'Interpreter', 'latex', 'FontSize', 20)
xlabel('Number of Particles ($N$)', 'Interpreter', 'latex', 'FontSize', 18)
ylabel('Run Time (s)', 'Interpreter', 'latex', 'FontSize', 18)
legend({'SPaM', 'PLSS'}, 'Interpreter', 'latex', 'FontSize', 16, 'Location', 'best')

% Error plot (linear)
figure(3)
clf
plot(Np_list, DDCErr(1, :))
hold on
plot(Np_list, DDCErr(2, :))
title('RMSE (linear)', 'Interpreter', 'latex', 'FontSize', 20)
xlabel('Number of Particles ($N$)', 'Interpreter', 'latex', 'FontSize', 18)
ylabel('RMSE', 'Interpreter', 'latex', 'FontSize', 18)
legend({'SPaM', 'PLSS'}, 'Interpreter', 'latex', 'FontSize', 16, 'Location', 'best')

% Error plot (loglog)
figure(4)
clf
loglog(Np_list, DDCErr(1, :))
hold on
loglog(Np_list, DDCErr(2, :))
title('RMSE (loglog)', 'Interpreter', 'latex', 'FontSize', 20)
xlabel('Number of Particles ($N$)', 'Interpreter', 'latex', 'FontSize', 18)
ylabel('RMSE', 'Interpreter', 'latex', 'FontSize', 18)
legend({'SPaM', 'PLSS'}, 'Interpreter', 'latex', 'FontSize', 16, 'Location', 'best')
