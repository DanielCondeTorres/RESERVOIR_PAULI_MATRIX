function truncate_operator(O::Operator, M::Int)
    # Si el número de términos ya es menor o igual a M, no hacemos nada
    if length(O) <= M
        return O
    end
    
    # 1. Convertimos el diccionario a un array de pares (PauliString => Coeficiente)
    # 2. Ordenamos por el valor absoluto del coeficiente de mayor a menor (rev=true)
    sorted_terms = sort(collect(O), by = x -> abs(x.second), rev=true)
    
    # Creamos un nuevo Operador solo con los M términos más grandes
    truncated_O = Operator()
    for i in 1:M
        truncated_O[sorted_terms[i].first] = sorted_terms[i].second
    end
    
    return truncated_O
end