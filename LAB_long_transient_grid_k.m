function [gridk,each_k_step] = LAB_long_transient_grid_k(Tfft,f0range,f0dotrange,nb)

% This function computes the grid on the coefficient k of the power law
% fdot \propto f^n. The grid is computed using an analytical formula. 
%
% Input arguments:
%
% Tfft: fft length
% f0range: frequency range to be explored [Hz]
% f0dotrange: spin-down range [Hz/s]
% tobs: total observation time [s]. Needed for tests (not for the grid
% computation)

% Output argument:
%
% gridk: grid on the coefficient k
%
% Use example: gridk = long_transient_grid_k(1,[500, 2000],[-1e-4 -1],1e5); 
% C. Palomba and the Rome CW group (Ocotber 2017)

% min and max search frequency
f0min=min(f0range); 
f0max=max(f0range);
% min and max search spin-down 
f0dot_min=min(abs(f0dotrange));
f0dot_max=max(abs(f0dotrange));

% Here we make a loop on the braking index. In fact the loop should go 
% over the braking index grid computed by long_transient_grid_n--> to be implemented 
nn=1;
for i=1:nn % loop on braking index 
    %n(i)=2+(i-1)*5/(nn-1); % current braking index value
    %nb=n(i);
    kmin=f0dot_min/f0max^nb; % corresponding kmin
    kmax=f0dot_max/f0min^nb; % corresponding kmax
    fprintf('k range: %e %e\n',kmin,kmax);
    newk=kmax; % start from the maximum k
    j=1;
    clear k
    while newk >= kmin % loop as long as the k goes below kmin
        dk=newk*((1+1./(Tfft*f0max))^(-nb)-1); % step in k (note: it's function of k. 
        each_k_step(j)=dk;
        % This relation can be computed analytically.)
        % k(j) contains the values of k for the current braking index.
        if j==1
            k(j)=newk;
        else
            k(j)=newk+dk; % next k
        end
        newk=k(j);
        j=j+1;
    end
    nk(i)=length(k); % number of k values in the range kmin-kmax for the current braking index
    gridk=k; % grid on k for the i-th braking index value
    
    %%%%%%%%%%%%%%%%
    % Test on the k range. We take random values of f0 and f0dot in the
    % search range given in input and compute the final frequencies for the
    % minumum and maximum k values. Final frequencies too low mean that
    % probably the range of k is not completely meaningful and should be
    % reduced. To be completed.
%     randfdot=10.^log10randfdot;
%     ftobsrand=randf./(1.+(nb-1)*randfdot./randf.*tobs).^(1./(nb-1));
%     kkk=find(ftobsrand>=20.0);
%     randfdot_sel=randfdot(kkk);
%     randf_sel=randf(kkk);
%     ftobsrand_sel=ftobsrand(kkk);
%     %figure;plot_triplets_noauto(randfdot_sel,randf_sel,ftobsrand_sel,'.','hot');
%     ftobskmin(i)=f0max/(1+(nb-1)*kmin*f0max^(nb-1)*tobs)^(1/(nb-1)); % maximum frequency at tobs
%     ftobskmax(i)=f0min/(1+(nb-1)*kmax*f0min^(nb-1)*tobs)^(1/(nb-1)); % minimum frequency at tobs
    %%%%%%%%%%%%%%%%%% 
end
gridk=flip(gridk); %smallest first
each_k_step=abs(flip(each_k_step)); %so corresponding dk is to the smallest gridk, etc.

fprintf('Total number of k values %d\n',length(gridk));
