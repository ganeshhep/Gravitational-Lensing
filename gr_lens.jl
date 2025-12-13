using Plots

function lensing_plots(D_s, D_d, D_ds, R_s;
                       eta_tilde,
                       nphi)

    # Lensing scales

    alpha_0 = sqrt((2 * R_s * D_ds) / (D_d * D_s))
    eta_0 = alpha_0 * D_s
    xi_0 = alpha_0 * D_d

    # Source distribution (circle)
 
    phi = range(0, 2pi, length = nphi)

    x_src = @. eta_tilde * eta_0 * cos(phi)
    y_src = @. eta_tilde * eta_0 * sin(phi)

    # Image positions
 
    xi_plus  = 0.5 * (eta_tilde + sqrt(eta_tilde^2 + 4)) * xi_0
    xi_minus = 0.5 * (eta_tilde - sqrt(eta_tilde^2 + 4)) * xi_0

    x_img1 = @. xi_plus * cos(phi)
    y_img1 = @. xi_plus * sin(phi)

    x_img2 = @. xi_minus * cos(phi)
    y_img2 = @. xi_minus * sin(phi)

    # Source Plot in the Source Plane

    p = scatter(
        x_src, y_src,
        markersize = 2,
        label = "Source",
        xlabel = "x",
        ylabel = "y",
        title = "Source Plane",
        size = (2000, 1500),
        dpi = 300
    )

    savefig(p, "source.png")

    # Image Plot in the Lens Plane
 
    q = scatter(
        x_img1, y_img1,
        markersize = 2,
        label = "Outer Image",
        xlabel = "x",
        ylabel = "y",
        title = "Lens Plane",
        size = (2000, 1500),
        dpi = 300
    )

    scatter!(
        q,
        x_img2, y_img2,
        markersize = 2,
        label = "Inner Image"
    )

    savefig(q, "image.png")
    
end

# Distances in Parsec

D_s = 2e3
D_d = 1e3
D_ds = 1e3
R_s = 5.4776e-12

lensing_plots(D_s, D_d, D_ds, R_s;
              eta_tilde = 3.0,
              nphi = 1500)


