thresh=3.5;
l=210;

[eff,eff_dff,eff_gen]=deal(zeros(1,25));
for i=1:25
    eff(i)=sum(crmat(i,:)>thresh)/l;
    eff_gen(i)=sum(crmat_gen(i,:)>thresh)/l;
    eff_dff(i)=sum(crmat_dff(i,:)>thresh)/l;
end

figure,semilogx(h0, eff,'.-','MarkerSize',6)
hold on
semilogx(h0, eff_dff,'.-','MarkerSize',6)
semilogx(h0, eff_gen,'.-','MarkerSize',6)
title('Efficiency with threshold=3.5')
xlabel('h0')
ylabel('Efficiency')
ylim([-0.07 1.07])
xlim([h0(1) 4E-23])
grid on
legend('normal','custom filter','generalised filter','Location','northwest')
hold off