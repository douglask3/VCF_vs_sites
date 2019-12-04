is_p_star <- function(P) {
    if (P > 0.1) out = ' '
    else if (P > 0.05) out = '.'
    else if (P > 0.01) out = '*'
    else if (P > 0.001) out = '**'
    else out = '***'
    out
}