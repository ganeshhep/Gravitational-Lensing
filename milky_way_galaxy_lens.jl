using Plots
plotlyjs()
using Contour
using FiniteDifferences

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

function ψ_bulge(m_b, x_1, x_2)
    x = sqrt(x_1^2 + x_2^2)

    ψ₁₁ = scaled_deflection_prime(m_b, x) * (x_1^2/x^2) + scaled_deflection(m_b, x) * (x_2^2/x^3)
    ψ₂₂ = scaled_deflection_prime(m_b, x) * (x_2^2/x^2) + scaled_deflection(m_b, x) * (x_1^2/x^3)
    ψ₁₂ = scaled_deflection_prime(m_b, x) * ((x_1 * x_2)/x^2) - scaled_deflection(m_b, x) * ((x_1 * x_2)/x^3)
    
    return ψ₁₁, ψ₂₂, ψ₁₂
end

function ψ_disc(m_d, a, b, x_1, x_2)
    term_1 = (a + sqrt(x_2^2 + b^2))^2 - x_1^2
    term_2 = ((a + sqrt(x_2^2 + b^2))^2 + x_1^2)^2
    ψ₁₁ = m_d * (term_1/term_2)

    term_3 = -2 * x_1 * x_2 * (a + sqrt(x_2^2 + b^2))
    term_4 = sqrt(x_2^2 + b^2) * (x_1^2 + (a + sqrt(x_2^2 + b^2))^2)^2
    ψ₁₂ = m_d * (term_3/term_4)

    term_5 = a * sqrt(x_2^2 + b^2) + 2 * x_2^2 + b^2
    term_6 = (x_1^2 + (a + sqrt(x_2^2 + b^2))^2) * (x_2^2 + b^2)
    term_7 = x_2^2 * (a + sqrt(x_2^2 + b^2)) * (4 * a * sqrt(x_2^2 + b^2) + 3 * (x_2^2 + b^2) + (x_1^2 + a^2))
    term_8 = (x_1^2 + (a + sqrt(x_2^2 + b^2))^2)^2 * (x_2^2 + b^2)^1.5
    ψ₂₂ = m_d * (term_5/term_6) - m_d * (term_7/term_8)
    
    return ψ₁₁, ψ₂₂, ψ₁₂
end

function ψ_halo(ρ_h, λ, x_1, x_2)
    x = sqrt(x_1^2 + x_2^2)

    term_1 = sqrt(1 + λ^2 * x^2) + 1 + λ^2 * x_2^2
    term_2 = sqrt(1 + λ^2 * x^2) * (1 + sqrt(1 + λ^2 * x^2))^2
    ψ₁₁ = 2 * ρ_h * (term_1/term_2)

    term_3 = sqrt(1 + λ^2 * x^2) + 1 + λ^2 * x_1^2
    ψ₂₂ = 2 * ρ_h * (term_3/term_2)

    term_4 = - λ^2 * x_1 * x_2
    ψ₁₂ = 2 * ρ_h * (term_4/term_2)
    
    return ψ₁₁, ψ₂₂, ψ₁₂
end

function ψ_milky_way(m_b, m_d, a, b, ρ_h, λ, x_1, x_2)
    ψ_b_11, ψ_b_22, ψ_b_12 = ψ_bulge(m_b, x_1, x_2)
    ψ_d_11, ψ_d_22, ψ_d_12 = ψ_disc(m_d, a, b, x_1, x_2)
    ψ_h_11, ψ_h_22, ψ_h_12 = ψ_halo(ρ_h, λ, x_1, x_2)
    
    ψ₁₁ = ψ_b_11 + ψ_d_11 + ψ_h_11
    ψ₂₂ = ψ_b_22 + ψ_d_22 + ψ_h_22
    ψ₁₂ = ψ_b_12 + ψ_d_12 + ψ_h_12

    return ψ₁₁, ψ₂₂, ψ₁₂
end

function α_bulge(m_b, x_1, x_2)
    x = sqrt(x_1^2 + x_2^2)

    α₁ = scaled_deflection(m_b, x) * (x_1/x)
    α₂ = scaled_deflection(m_b, x) * (x_2/x)

    return α₁, α₂
end

function α_disc(m_d, a, b, x_1, x_2)
    term_1 = x_1
    term_2 = x_1^2 + (a + sqrt(x_2^2 + b^2))^2
    α₁ = m_d * (term_1/term_2)

    term_3 = x_2 * (a + sqrt(x_2^2 + b^2))
    term_4 = sqrt(x_2^2 + b^2) * (x_1^2 + (a + sqrt(x_2^2 + b^2))^2)
    α₂ = m_d * (term_3/term_4)

    return α₁, α₂
end

function α_halo(ρ_h, λ, x_1, x_2)
    x = sqrt(x_1^2 + x_2^2)

    term = 1 + sqrt(1 + λ^2 * x^2)

    α₁ = 2 * ρ_h * (x_1/term)
    α₂ = 2 * ρ_h * (x_2/term)

    return α₁, α₂
end

function α_milky_way(m_b, m_d, a, b, ρ_h, λ, x_1, x_2)
    α_b_1, α_b_2 = α_bulge(m_b, x_1, x_2)
    α_d_1, α_d_2 = α_disc(m_d, a, b, x_1, x_2)
    α_h_1, α_h_2 = α_halo(ρ_h, λ, x_1, x_2)

    α₁ = α_b_1 + α_d_1 + α_h_1
    α₂ = α_b_2 + α_d_2 + α_h_2

    return α₁, α₂
end

function det(m_b, m_d, a, b, ρ_h, λ, x_1, x_2)
    ψ₁₁, ψ₂₂, ψ₁₂ = ψ_milky_way(m_b, m_d, a, b, ρ_h, λ, x_1, x_2)

    κ = 0.5 * (ψ₁₁ + ψ₂₂)
    γ₁ = 0.5 * (ψ₁₁ - ψ₂₂)
    γ₂ = ψ₁₂
    γ = sqrt(γ₁^2 + γ₂^2)
    
    return (1 - κ)^2 - γ^2
end

function psrc_img(m_b, m_d, a, b, ρ_h, λ, y1_src, y2_src, min, max, N, tol)
    x1_vals = x2_vals = range(min, max, N)

    x1_sol = Float64[]
    x2_sol = Float64[]

    for x1 in x1_vals, x2 in x2_vals
        α₁ = α_milky_way(m_b, m_d, a, b, ρ_h, λ, x1, x2)[1]
        α₂ = α_milky_way(m_b, m_d, a, b, ρ_h, λ, x1, x2)[2]

        y1 = x1 - α₁
        y2 = x2 - α₂

        dist = sqrt((y1 - y1_src)^2 + (y2 - y2_src)^2)

        if dist < tol
           push!(x1_sol, x1)
           push!(x2_sol, x2)
        end
    end

    return x1_sol, x2_sol
end

function mwg_lens(m_b, m_d, a, b, ρ_h, λ, y1_src, y2_src, min, max, N, tol)
    x1_val = x2_val = range(min, max, N)

    det_vals = [det(m_b, m_d, a, b, ρ_h, λ, x1_c, x2_c) for x1_c in x1_val, x2_c in x2_val]

    c = Contour.contour(x1_val, x2_val, det_vals, 0.0)

    cd = central_fdm(5, 1)

    ψ₁₁_num = [cd(x_1 -> α_milky_way(m_b, m_d, a, b, ρ_h, λ, x_1, x_2)[1], x_1) for x_1 in x1_val, x_2 in x2_val]
    ψ₂₂_num = [cd(x_2 -> α_milky_way(m_b, m_d, a, b, ρ_h, λ, x_1, x_2)[2], x_2) for x_1 in x1_val, x_2 in x2_val]
    ψ₁₂_num = [cd(x_1 -> α_milky_way(m_b, m_d, a, b, ρ_h, λ, x_1, x_2)[2], x_1) for x_1 in x1_val, x_2 in x2_val]

    κ = @. 0.5 * (ψ₁₁_num + ψ₂₂_num)
    γ₁ = @. 0.5 * (ψ₁₁_num - ψ₂₂_num)
    γ₂ = ψ₁₂_num
    γ = @. sqrt(γ₁^2 + γ₂^2)

    det_vals_num = @. (1 - κ)^2 - γ^2

    # println(det_vals[1][1])
    # println(det_vals_num[1][1])

    cnum = Contour.contour(x1_val, x2_val, det_vals_num, 0.0)

    plt1 = plot(; legend = true)
    plt2 = plot(; legend = true)
    plt3 = plot(; legend = true)
    
    # Analytical Solutions
    # for line in lines(c)
    #     x1, x2 = coordinates(line)

    #     l = length(x1)

    #     α₁ = Float64[]
    #     α₂ = Float64[]

    #     for i in 1:l
    #         ang1, ang2 = α_milky_way(m_b, m_d, a, b, ρ_h, λ, x1[i], x2[i])

    #         push!(α₁, ang1)
    #         push!(α₂, ang2)

    #     end

    #     y1 = @. x1 - α₁
    #     y2 = @. x2 - α₂

    #     scatter!(plt1, x1, x2, markersize = 1.0, label = "a = $a, b = $b, ρ = $ρ_h", xlabel = "x₁", ylabel = "x₂", title = "Milky Way Galaxy Lens : Critical Curves",
    #              size = (1000, 1000), dpi = 300)

    #     scatter!(plt2, y1, y2, markersize = 1.0, label = "a = $a, b = $b, ρ = $ρ_h", xlabel = "y₁", ylabel = "y₂", title = "Milky Way Galaxy Lens : Caustic Curves",
    #              size = (1000, 1000), dpi = 300)

    # end

    # Numerical Solutions
    for line_num in lines(cnum)
        x1_num, x2_num = coordinates(line_num)

        lnum = length(x1_num)

        α₁ = Float64[]
        α₂ = Float64[]

        for i in 1:lnum
            ang1, ang2 = α_milky_way(m_b, m_d, a, b, ρ_h, λ, x1_num[i], x2_num[i])

            push!(α₁, ang1)
            push!(α₂, ang2)

        end

        y1_num = @. x1_num - α₁
        y2_num = @. x2_num - α₂

        scatter!(plt1, x1_num, x2_num, markersize = 1.0, label = "a = $a, b = $b, ρ = $ρ_h", xlabel = "x₁", ylabel = "x₂", title = "Milky Way Galaxy Lens : Critical Curves",
                 size = (1000, 1000), dpi = 300, markercolor = "black")

        scatter!(plt2, y1_num, y2_num, markersize = 1.0, label = "a = $a, b = $b, ρ = $ρ_h", xlabel = "y₁", ylabel = "y₂", title = "Milky Way Galaxy Lens : Caustic Curves",
                 size = (1000, 1000), dpi = 300, markercolor = "black")

    end
    
    x1s, x2s = psrc_img(m_b, m_d, a, b, ρ_h, λ, y1_src, y2_src, min, max, N, tol)

    scatter!(plt3, x1s, x2s, label = "Image", xlabel = "x₁", ylabel = "x₂", title = "Milky Way Galaxy Lens : Point Source Images",
             size = (1000, 1000), dpi = 300, markercolor = "black")
    scatter!(plt3, [y1_src], [y2_src], label = "Source", markercolor = "blue")

    display(plt1)
    display(plt2)
    display(plt3)

end

a = 9.3
b = 0.37
m_b = 10.2
m_d = 30.0
ρ_h = 0.16
r_0 = 0.7
r_c = 6.0
λ = r_0/r_c
y1_src = 0.5
y2_src = 0.5
min = - 30.0
max = 30.0
N = 15000
tol = 1e-1

mwg_lens(m_b, m_d, a, b, ρ_h, λ, y1_src, y2_src, min, max, N, tol)
