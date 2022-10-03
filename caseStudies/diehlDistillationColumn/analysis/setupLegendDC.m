function setupLegendDC(ax)
% adds legend to axes
l = {'Interpreter','latex','location','southeast'};
legend(ax(1),{'$\overline{T}_{28,s}$','$\overline{T}_{14,s}$','$\overline{T}_{28,p}^{\ast}$','$\overline{T}_{14,p}^{\ast}$'},l{:})
legend(ax(2),{'$J_{p}$','$J_{p}^{\ast}$'},l{:})
legend(ax(3),{'$G_{1,p}$','$G_{2,p}$','$G_{1,p}^{\ast}$','$G_{2,p}^{\ast}$'},l{:})

end