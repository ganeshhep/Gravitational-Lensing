using Plots
# using MeshGrid
plotlyjs()
# default(show=true)
using Contour

function chi(x)
    if x > 1
        arg_1 = sqrt(x^2 - 1)
        return atan(arg_1)/arg_1
    elseif x < 1
        arg_2 = sqrt(1 - x^2)
        return atanh(arg_2)/arg_2
    else
        return 1.0
    end  
end

function chi_prime(x)
    term_1 = 1 - x^2 * chi(x)
    term_2 = x * (x^2 - 1)
    return term_1/term_2
end

function scaled_deflection(m_b, x)
    term_1 = x * (chi(x) - 1)
    term_2 = 1 - x^2
    return m_b * (term_1/term_2)
end

function scaled_deflection_prime(m_b, x)
    term_1 = ((1 + x^2)/(1 - x^2)^2) * (chi(x) - 1)
    term_2 = (x/(1 - x^2)) * chi_prime(x)
    
    return m_b * (term_1 + term_2)
end

function mean_surf_density(m_b, x) 
    return scaled_deflection(m_b, x)/x
end

function mean_surf_density_prime(m_b, x)
    term_1 = - scaled_deflection(m_b, x)/x^2 
    term_2 = scaled_deflection_prime(m_b, x)/x
    return term_1 + term_2 
end

function det(m_b, x_1, x_2, κₚ, γₚ)
    x = sqrt(x_1^2 + x_2^2)

    Γ₁ = κₚ + γₚ
    Γ₂ = κₚ - γₚ

    term_1 = 1 - mean_surf_density(m_b, x) - Γ₁
    term_2 = 1 - mean_surf_density(m_b, x) - Γ₂
    term_3 = (1 - mean_surf_density(m_b, x) - Γ₁) * (x_2^2/x) * mean_surf_density_prime(m_b, x) 
    term_4 = (1 - mean_surf_density(m_b, x) - Γ₂) * (x_1^2/x) * mean_surf_density_prime(m_b, x) 

    return term_1 * term_2 - term_3 - term_4
end

function critical_and_caustic_curves(m_b, min, max, N, κₚ, γₚ)
    Γ₁ = κₚ + γₚ
    Γ₂ = κₚ - γₚ
    
    x1_val = x2_val = range(min, max, N)

    det_vals = [det(m_b, x1_c, x2_c, κₚ, γₚ) for x1_c in x1_val, x2_c in x2_val]

    c = Contour.contour(x1_val, x2_val, det_vals, 0.0)

    plt1 = plot(; legend = true)
    plt2 = plot(; legend = true)

    for line in lines(c)
        x1, x2 = coordinates(line)

        x = @. sqrt(x1^2 + x2^2)

        y1 = @. (1 - mean_surf_density(m_b, x)) * x1 - Γ₁ * x1
        y2 = @. (1 - mean_surf_density(m_b, x)) * x2 - Γ₂ * x2

        scatter!(plt1, x1, x2, markersize = 1.0, label = "κₚ = $κₚ, γₚ = $γₚ", xlabel = "x₁", ylabel = "x₂", title = "Perturbed Hernquist Lens : Critical Curves",
                 size = (1000, 1000), dpi = 300)

        scatter!(plt2, y1, y2, markersize = 1.0, label = "κₚ = $κₚ, γₚ = $γₚ", xlabel = "y₁", ylabel = "y₂", title = "Perturbed Hernquist Lens : Caustic Curves",
                 size = (1000, 1000), dpi = 300)

    end

    display(plt1)
    display(plt2)
    
end

m_b = 1.0
min = -2.0
max = 2.0
N = 10000
κₚ = 0.0
γₚ = 0.2

critical_and_caustic_curves(m_b, min, max, N, κₚ, γₚ)