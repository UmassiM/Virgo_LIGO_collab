function out = Efficiency(crvec, thresh)

out=(sum(crvec>thresh)/length(crvec))*100;


end

