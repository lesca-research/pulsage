clear all
close all

file_name='tfa_sample_data.txt'; % select the file to process 

X=importdata(file_name); % get the data from the file

X=X.data; % get only the numerical values, not the column-header text
t=X(:,1); % first column is the time 
p1=X(:,2); % second column is the ABP
v1=X(:,3); % third column is CBFV left
v2=X(:,4); % fourth column is the CBFV right
fs=1/mean(diff(t)); % find the sampling frequency from the time-base

figure % plot the raw sigals
plot(t,p1,'k',t,v1,'k:');
title(file_name,'interpreter','none');
xlabel ('time (s)');
ylabel('mmHg, cm/s');
legend('ABP','CBFV');
title ('Raw signals');
axis tight

params.plot_title=[make_title(file_name), ',left CBFV'];% set the title of plots to the file-name 
%(make_title just makes sure that underscores do not become subscripts in the title). Leave all other parameters as default in tfa_car.m
tfa_out=tfa_car(p1,v1,fs,params); % apply tfa to the left CBFV channel (v1)
tfa_out

params.plot_title=[make_title(file_name), ',right CBFV'];% 
tfa_out=tfa_car(p1,v2,fs,params); % apply tfa to the right CBFV channel (v2)
tfa_out
