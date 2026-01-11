using Plots
plotlyjs()
using Contour

function det(m_d, a, b, x_1, x_2)
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

    return (1 - ψ₁₁) * (1 - ψ₂₂) - ψ₁₂^2
end

function α(m_d, a, b, x_1, x_2)
    term_1 = x_1
    term_2 = x_1^2 + (a + sqrt(x_2^2 + b^2))^2
    α₁ = m_d * (term_1/term_2)

    term_3 = x_2 * (a + sqrt(x_2^2 + b^2))
    term_4 = sqrt(x_2^2 + b^2) * (x_1^2 + (a + sqrt(x_2^2 + b^2))^2)
    α₂ = m_d * (term_3/term_4)

    return α₁, α₂
end

function critical_and_caustic_curves(m_d, a, b, min, max, N)
    x1_val = range(min, max, N)
    x2_val = range(min, max, N)
    
    x1_val = x2_val = range(min, max, N)

    det_vals = [det(m_d, a, b, x1_c, x2_c) for x1_c in x1_val, x2_c in x2_val]

    c = Contour.contour(x1_val, x2_val, det_vals, 0.0)

    plt1 = plot(; legend = true)
    plt2 = plot(; legend = true)

    for line in lines(c)
        x1, x2 = coordinates(line)

        l = length(x1)

        α₁ = Float64[]
        α₂ = Float64[]

        for i in 1:l
            ang1, ang2 = α(m_d, a, b, x1[i], x2[i]) 

            push!(α₁, ang1)
            push!(α₂, ang2)

        end

        y1 = @. x1 - α₁
        y2 = @. x2 - α₂

        scatter!(plt1, x1, x2, markersize = 1.0, label = "a = $a, b = $b", xlabel = "x₁", ylabel = "x₂", title = "Miyamoto Nagai Lens : Critical Curves",
                 size = (1000, 1000), dpi = 300)

        scatter!(plt2, y1, y2, markersize = 1.0, label = "a = $a, b = $b", xlabel = "y₁", ylabel = "y₂", title = "Miyamoto Nagai Lens : Caustic Curves",
                 size = (1000, 1000), dpi = 300)

    end

    display(plt1)
    display(plt2)

end

m_d = 1.0
a = 0.3
b = 0.5
min = -2.0
max = 2.0
N = 10000

critical_and_caustic_curves(m_d, a, b, min, max, N)
