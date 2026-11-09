close all
clear

names = ["u_stab1", "u_stab2", "u_unstab1", "u_unstab2"];
figs = ["sh_sol_stable", "sh_sol_2stable", "sh_sol_unstable", "sh_sol_2unstable"];
ds = [7, 12, 9, 15]; % the domain is [-d,d]^2
views = [3, 2];

for k = 1:4
    % each line of the file: n1 n2 real(a_n) imag(a_n)
    data = load(names(k) + ".txt");
    n1 = data(:,1);
    n2 = data(:,2);
    a = data(:,3) + 1i*data(:,4);
    d = ds(k);
    x = -d:0.1:d;
    [X1, X2] = meshgrid(x);
    % u_0(x) = sum_n a_n exp(i*pi/d*(n1*x1 + n2*x2))
    U = zeros(size(X1));
    for j = 1:length(a)
        U = U + a(j)*exp(1i*pi/d*(n1(j)*X1 + n2(j)*X2));
    end
    U = real(U);

    for v = 1:2
        fig = figure;
        surf(x, x, U)
        shading interp
        colormap(jet);
        axis tight

        set(gca,'FontSize',12)
        xlabel('$$x_1$$', 'Interpreter', 'latex', 'FontSize', 20)
        ylabel('$$x_2$$', 'Interpreter', 'latex', 'FontSize', 20)
        zlabel('$$u_0(x)$$', 'Interpreter', 'latex', 'FontSize', 15)
        colorbar
        view(views(v))
        axis square
        xticks([-d 0 d])
        yticks([-d 0 d])

        position = fig.Position;
        position(4) = min(position(3:4));
        position(3) = round(1.2*position(4));
        fig.Position = position;
        set(gca, 'Position', [0.17 0.165 0.617 0.76])
        theme(fig, "light")
        exportgraphics(fig, "figures/" + figs(k) + "_2D_view" + v + ".pdf", "ContentType", "image", "Resolution", 600, "BackgroundColor", "white")
        close(fig)
    end
end
