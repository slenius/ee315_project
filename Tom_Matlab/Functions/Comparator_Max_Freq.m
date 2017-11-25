function [ comp_dec_time, reg_time, recharge_time ] = Comparator_Max_Freq( vop, vom, vclock, num_clock_sims )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for i=1:num_clock_sims
    clearvars -EXCEPT comp_dec_time reg_time recharge_time vop vom vclock num_clock_sims i
    % first find total decision time. look for time it takes to reach
    % 0.9VDD after clock goes high
    vdiff.V=abs(vop.V(:,i)-vom.V(:,i));
    vdiff.time=vop.time(:,i);
    search_1=0.5;
    search_2=0.9;
    Vdd=1.2;
    % simulation set up so clock pulse begins at 10ns. Total clock high
    % time should last no longer than 5ns
    time_start=10e-9;
    time_stop=14.9e-9;
    % find index of clock high and middle
    [~,index_clock_start]=min(abs(vclock.time-time_start));
    [~,index_clock_middle]=min(abs(vclock.time-time_stop));
    % find time to 0.9 VDD in the differential output
    diff_stop_val=0.9*Vdd;
    [~,index_diff_high]=min(abs(vdiff.V(index_clock_start:index_clock_middle)-diff_stop_val));
    index_diff_high=index_diff_high+index_clock_start;
    comp_dec_time(i)=vclock.time(index_diff_high)-vclock.time(index_clock_start);
    
    % now find regeneration time. This is the time when the differential
    % output begins rising from 0.1VDD to 0.9VDD
    [~,index_diff_reg_start]=min(abs(vdiff.V(index_clock_start:index_clock_middle)-Vdd*0.1));
    index_diff_reg_start=index_diff_reg_start+index_clock_start;
    reg_time(i)=vclock.time(index_diff_high)-vclock.time(index_diff_reg_start);
    
    % now find time to recharge after clock goes low. This will be from
    % 15ns to the end of the dataset as defined in simulation
    time_clock_fall=15e-9;
    % find index of clock fall
    [~,index_clock_fall]=min(abs(vclock.time-time_clock_fall));
    % find time for comparator output to fall back to 0.1VDD
    [~,index_comp_reset]=min(abs(vdiff.V(index_clock_fall:end)-Vdd*0.1));
    index_comp_reset=index_comp_reset+index_clock_fall;
    
    recharge_time(i)=vclock.time(index_comp_reset)-vclock.time(index_clock_fall);
   
end


end

