function [peakss,indss,times] = LAB_create_peakmap_from_wnimg(img_in,thresh,peak_flag,zero_flag)



if isstruct(img_in)
    img=img_in.Z;
    t0=img_in.inix;
    f0=img_in.iniy;
    dt=img_in.dx;
    df=img_in.dy;
else
    img=img_in;
    t0=0;
    f0=0;
    dt=1;
    df=1;
end

[Nt,Nf]=size(img);
peakss=[];
indss=[];

for i_t=1:Nt
    ratio=img(i_t,:);  %IMP the sqrt(sqrt(2)) for simulated noise

    peaks_index=1; %index for peaks
    for i_f=2:(Nf-1)%1:length(ratio)     
        if peak_flag==1
            if ratio(i_f)>thresh %for all ratios calculated, check if it's first bigger than threshold %%% NEED THIS SETUP FOR THE HOUGH
                 if ratio(i_f)>ratio(i_f+1)
                     if ratio(i_f)>ratio(i_f-1) %check in middle of FFT, if bigger than both the one before and one after
                        peaks(1,peaks_index)=t0+dt*(i_t-1);
                        peaks(2,peaks_index)=f0+df*(i_f-1); %adding frequency from beginning
                        peaks(3,peaks_index)=ratio(i_f);
                        peaks_index=peaks_index+1;
                     end
                 end
            end            
        else
            peaks(1,peaks_index)=t0+dt*(i_t-1);
            peaks(2,peaks_index)=f0+df*(i_f-1); %adding frequency from beginning

            if zero_flag==0 %put zeros where the CR is not a local max, applied threshold, so that the size of each peakmap is the same

                if ratio(i_f)>thresh %for all ratios calculated, check if it's first bigger than threshold %%% NEED THIS SETUP FOR THE HOUGH
                    if ratio(i_f)>ratio(i_f+1)
                        if ratio(i_f)>ratio(i_f-1) %check in middle of FFT, if bigger than both the one before and one after
                            peaks(3,peaks_index)=ratio(i_f);
                        else
                            peaks(3,peaks_index)=0;
                        end
                    else
                        peaks(3,peaks_index)=0;
                    end
                else
                    peaks(3,peaks_index)=0;
                end

            else
                peaks(3,peaks_index)=ratio(i_f); %records CR if it's less than thresh, doesn't make it 0
            end
            peaks_index=peaks_index+1;
        end
            %if none of these are true, i = i+1, the loop advances
    end
    if ~exist('peakss','var')
       peaks=[];  
    end
    peakss=[peakss peaks];
    indss=[indss peaks_index];
end


times=t0+(0:(Nt-1))*dt;

end




