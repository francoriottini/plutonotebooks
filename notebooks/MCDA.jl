### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 8d16ed10-40f1-11ee-22c7-0d2050858a0c
using Pluto, PlutoUI, DataFrames, LinearAlgebra, StatsPlots, Plots, PlotlyJS, FileIO, Cairo, Compose, Images, HTTP, InteractiveUtils, Interact, WebIO

# ╔═╡ e5ce8ca5-ac4c-437f-b6e7-0ad8df985358
begin 
	proyectos = [
    "Corredor pacífico de integración logística Mesoamericana",
    "Conglomerado portuario del Caribe",
    "Integración eléctrica latinoamericana",
    "Integración transfronteriza andina",
    "Corredor bioceánico de capricornio",
    "Iniciativa de la Cuenca del Plata",
    "Corredor de conectividad MERCOSUR-Chile",
    "Conectividad bioceánica en la Patagonia",
    "Infraestructura digital y de datos para la integración social y productiva en ALC"]
	
	md""" Los proyectos a tratar son:
	
	-- Corredor pacífico de integración logística Mesoamericana

	-- Conglomerado portuario del Caribe

	-- Integración eléctrica latinoamericana

	-- Integración transfronteriza andina

	-- Corredor bioceánico de capricornio

	-- Iniciativa de la Cuenca del Plata

	-- Corredor de conectividad MERCOSUR-Chile

	-- Conectividad bioceánica en la Patagonia

	-- Infraestructura digital y de datos para la integración social y productiva en ALC"""
end

# ╔═╡ 8aa54581-5d3b-42e4-aaf8-8d80d00b14d2
begin 
	criterios = [
    "Impacto Económico", "Impacto Social", "Impacto Ambiental",
    "Cooperación", "Prioridad Política", "Madurez Operativa",
    "Ventaja Comparativa del BID", "Movilización de Recursos"]

	md""" Y los criterios:

	-- Impacto Económico, Impacto Social e Impacto Ambiental.

	-- Cooperación.

	-- Prioridad Política y Madurez Operativa.

	-- Ventaja Comparativa del BID y Movilización de Recursos."""

end
	

# ╔═╡ f6cb3e82-5f54-4e3c-a13d-50eeac715463
begin 
	criterios_principales = ["Impacto", "Cooperación", "Viabilidad", "Soporte Institucional"]
	md""""""

end

# ╔═╡ ff82c5c0-5139-4a4f-977c-f9193f639051
md"""### Armamos la tabla inicial de comparaciones entre criterios. Cada uno representa la comparción de i con j. Por ejemplo, "Impacto" es 3 veces más importante que "Cooperación"."""

# ╔═╡ 795eecd5-8784-4dff-b753-36662af88dd3
begin 
	subcriterios_impacto = ["Económico", "Social", "Ambiental"]
	md""""""

end

# ╔═╡ 9a55981b-81e4-4601-89e5-0dbb41ce8f36
begin 
	subcriterios_viabilidad = ["Prioridad Política", "Madurez Operativa"]
	md""""""

end

# ╔═╡ 5a14d156-41be-47d8-84d7-e89aae7142a8
md"""### Para poder calcular los autovalores y autovectores la matriz debe ser cuadrada y recíproca. Por eso antes de que sigamos tenemos que chequear que la matriz (las matrices) que armamos lo sean. Eso hace el paso siguiente."""

# ╔═╡ 9c116dd2-5e11-406c-83e4-1f3a7dd6d810
function check_matrix(matriz::Matrix{Float64})
    # Comprobar si la matriz es cuadrada
    if size(matriz, 1) != size(matriz, 2)
        return false, "La matriz no es cuadrada."
    end

    # Comprobar si la matriz es recíproca
    for i in 1:size(matriz, 1)
        for j in i+1:size(matriz, 2)
            if matriz[i, j] * matriz[j, i] != 1.0
                return false, "La matriz no es recíproca."
            end
        end
    end

    return true, "La matriz es cuadrada y recíproca."
end


# ╔═╡ 1e8d91bd-ae61-4bcb-ab94-f262ebce1319
function calcular_ponderadores(matriz::Matrix{Float64})
    # Suma de cada columna
    suma_columnas = sum(matriz, dims=1)

    # Dividir cada elemento de la matriz por la suma de su columna
    matriz_normalizada = matriz ./ suma_columnas

    # Sumar cada fila y dividir por el número de criterios para obtener los ponderadores
    ponderadores = sum(matriz_normalizada, dims=2) ./ size(matriz, 1)

    return ponderadores
end

# ╔═╡ 6922f54f-4b48-4e0d-97c1-8469c4550f5c
md"""# Además, chequeamos la consistencia siguiendo el metodo descripto en Pujadas et al. (2017). Todos los CR deben dar menores 0.1"""

# ╔═╡ acec33c0-9012-42f2-91e7-f18dbb5de7e4
function calcular_consistencia(matriz::Matrix{Float64})
    n = size(matriz, 1) # Tamaño de la matriz

    # Valores y vectores propios
    eigvals, _ = eigen(matriz)

    # Valor propio máximo
    λ_max = maximum(real.(eigvals))

    # Índice de Consistencia
    CI = (λ_max - n) / (n - 1)

    # Índice Aleatorio (según la tabla proporcionada)
    RI = [0.00, 0.00, 0.58, 0.9, 1.12, 1.24, 1.32, 1.41, 1.45, 1.51][n]

    # Razón de Consistencia
    CR = CI / RI

    return CR
end

# ╔═╡ 8c007656-930c-4410-8d7c-7815d43b64c7


# ╔═╡ 6137955d-0201-4e57-a6e7-7dbb70b90645
md"""El método de Saaty para calcular los pesos utilizando la matriz de eigenvectores y eigenvalores puede ser complejo y requerir un conocimiento más avanzado de álgebra matricial. Además, en este contexto de análisis MCDA, el método más sencillo de cálculo mediante la media geométrica de las filas es preferible debido a su simplicidad y facilidad de implementación. Aunque los resultados entre ambos métodos pueden ser muy similares, la opción más sencilla es más accesible y menos propensa a errores para aquellos que no estén familiarizados con conceptos matriciales más avanzados."""

# ╔═╡ f716db17-1bdd-4d4e-9fe2-a578613bfb03
md""" ## Calculamos los pesos para las matrices"""

# ╔═╡ ae3aa23c-0895-42ee-9989-7ef66459cdd1
function calcular_pesos(matriz::Matrix{Float64})
    n = size(matriz, 1)
    
    # Cálculo de la media geométrica de cada fila
    geom_means = prod(matriz, dims=2).^(1/n)
    
    # Suma de las medias geométricas
    total_geom_mean = sum(geom_means)
    
    # Normalización de los pesos
    pesos = geom_means / total_geom_mean
    
    return pesos
end

# ╔═╡ c938c87e-b87e-49db-a585-e3c6d80196f9
md"""# Criterios para el proyecto "Corredor pacífico de integración logística Mesoamericana"
"""

# ╔═╡ dfdf0b56-f9ad-427c-a600-5f56ce8596a0
md"""Armamos un vector con todos los valores asignados a cada criterio"""

# ╔═╡ b0c81f63-8e86-4fa9-a764-5d9f7b772ea6
md"""Calculamos la suma ponderada para las secciones con subcriterios"""

# ╔═╡ 2aab428f-f7ef-4a52-b928-13888e517db7
md""" Creamos un nuevo vector que tiene en cuenta las sumas ponderadas de las secciones con subcriterios"""

# ╔═╡ 3b9f9560-fc6e-488a-a7fd-dcb71710f2e4
md""" Ahora calculamos la suma ponderada final utilizando los pesos principales"""

# ╔═╡ a9e1f4a6-782b-44c2-9523-c9892731ccdf
md"""# Criterios para el proyecto "Conglomerado portuario del Caribe"
"""

# ╔═╡ c6e082bb-980d-4bb5-9e94-3d3197d6f699
md"""Armamos un vector con todos los valores asignados a cada criterio"""

# ╔═╡ 62481b55-70f2-4525-811e-3d88a8bd41e1
md"""Calculamos la suma ponderada para las secciones con subcriterios"""

# ╔═╡ 6043c869-030e-4728-bfe0-731498a8c9ad
md"""# Criterios para el proyecto "Integración eléctrica latinoamericana"
"""

# ╔═╡ 08b56360-6a2f-4ee3-ab91-ad4d7c74a18b
md"""Calculamos la suma ponderada para las secciones con subcriterios"""

# ╔═╡ 7adc6966-4a50-414a-bbb5-38efc9ffd51c
md"""# Criterios para el proyecto "Integración transfronteriza andina"
"""

# ╔═╡ fa088abb-5d03-472e-b693-4eca3efc58b8
md"""Armamos un vector con todos los valores asignados a cada criterio"""

# ╔═╡ 93fa0b3e-3592-4140-aea3-237a414f9eaa
md"""Calculamos la suma ponderada para las secciones con subcriterios"""

# ╔═╡ 8c598b3c-f564-4d8f-b007-993f520b1589
md"""# Criterios para el proyecto "Corredor bioceánico de capricornio"
"""

# ╔═╡ f75ffe35-159c-4e9e-8956-16e55816c094
md"""Calculamos la suma ponderada para las secciones con subcriterios"""

# ╔═╡ 7e489c8d-e253-40ed-ba64-fb0672f5b3bb
md"""# Criterios para el proyecto "Iniciativa de la Cuenca del Plata"
"""

# ╔═╡ 0e6ad1ac-a123-42dc-865a-d664ced846f0
md"""Calculamos la suma ponderada para las secciones con subcriterios"""

# ╔═╡ 9f2cc85b-b55e-4be5-8d68-b442a20b150d
md"""# Criterios para el proyecto "Corredor de conectividad MERCOSUR-Chile"
"""

# ╔═╡ 27aa5269-df7d-4d84-b671-32167ca254a8
proyecto7 = proyectos[7]

# ╔═╡ 355aebc0-858d-4953-8a49-3e0315954e92
md"""Armamos un vector con todos los valores asignados a cada criterio"""

# ╔═╡ ffd6ca37-0c54-49d7-8e0d-9f3f0601f0bf
md"""Calculamos la suma ponderada para las secciones con subcriterios"""

# ╔═╡ e863b703-5476-4297-87bb-aa378caf8070
md"""# Criterios para el proyecto "Conectividad bioceánica en la Patagonia"
"""

# ╔═╡ 95768a04-cc66-422b-9a84-c67da26ed7cc
import PlutoUI: combine

# ╔═╡ bf7490a5-400c-4885-9be3-323320e6c677
function princip_pond_select(criterios::Dict)

	return combine() do Child
		
		inputs = [
			md""" $(v): $(
				Child(k, NumberField(0:0.05:100, default=0))
			)"""
			
			for (k, v) in criterios
		]
		
		md"""
		#### Comparación de criterios principales: Indique cuan más importante es el primer criterio en comparación con el segundo.
		$(inputs)
		"""
	end
end


# ╔═╡ 784194ab-5637-4458-9f6c-4d182c25ed41
@bind valor princip_pond_select(Dict(
    :I_vs_E => "Impacto versus Cooperacion", 
    :I_vs_V => "Impacto versus Viabilidad", 
    :I_vs_S => "Impacto versus Soporte", 
    :C_vs_V => "Cooperacion versus Viabilidad", 
    :C_vs_S => "Cooperacion versus Soporte", 
    :V_vs_S => "Viabilidad versus Soporte"
))

# ╔═╡ d1a1cbbf-76ce-46bd-8a8f-f0dfb4ce996b
begin 
	c1 = valor.I_vs_E
	c2 = valor.I_vs_V
	c3 = valor.I_vs_S
	c4 = valor.C_vs_V
	c5 = valor.C_vs_S
	c6 = valor.V_vs_S
	md""""""

end

# ╔═╡ 29beb7d0-04f0-4620-b5c3-644583f066e4
begin 
	matriz_criterios_principales = [
    1   c1   c2   c3;
    1/c1 1   c4 c5;
    1/c2 1/c4   1   c6;
    1/c3 1/c5 1/c6 1
]
	md""""""

end

# ╔═╡ 81c4e439-7cf5-4874-9ce6-33f3fdd66515
begin
    tabla_criterios = DataFrame(matriz_criterios_principales, Symbol.(criterios_principales))
	# Añadimos una columna con los nombres de los criterios
    tabla_criterios[!, "Comparación"] = criterios_principales

    # Reordenamos las columnas para que "Comparación" sea la primera
    tabla_criterios = tabla_criterios[!, ["Comparación"; criterios_principales]]
end

# ╔═╡ 281309e0-e252-4e06-bd49-c216e3ad0502
function impacto_pond_select(criterios::Dict)

	return combine() do Child
		
		inputs = [
			md""" $(v): $(
				Child(k, NumberField(0:0.05:100, default=0))
			)"""
			
			for (k, v) in criterios
		]
		
		md"""
		#### Comparación de subcriterios de Impacto: Indique cuan más importante es el primer subcriterio en comparación con el segundo.
		$(inputs)
		"""
	end
end

# ╔═╡ 530978c0-5aae-46cd-b7c3-987d17442311
@bind valor_impacto impacto_pond_select(Dict(
    :E_vs_S => "Impacto Económico versus Impacto Social", 
    :E_vs_A => "Impacto Económico versus Impacto Ambiental", 
    :S_vs_A => "Impacto Social versus Impacto Ambiental"
))

# ╔═╡ ef636cd9-1d27-43fc-9323-cadc6a4df699
begin 
	sci1 = valor_impacto.E_vs_S
	sci2 = valor_impacto.E_vs_A
	sci3 = valor_impacto.S_vs_A
	md""""""

end

# ╔═╡ 9de4a365-b0f1-46d4-991c-42188e4babcc
begin 
	matriz_subcriterios_impacto = [
    1   sci1   sci2;
    1/sci1 1   sci3;
    1/sci2 1/sci3 1
]
	md""""""

end

# ╔═╡ eb74d652-7fe1-4cd3-9013-a11d3efb7f49
begin
    tabla_subcriterios_impacto = DataFrame(matriz_subcriterios_impacto, Symbol.(subcriterios_impacto))
	# Añadimos una columna con los nombres de los criterios
    tabla_subcriterios_impacto[!, "Comparación"] = subcriterios_impacto

    # Reordenamos las columnas para que "Comparación" sea la primera
    tabla_subcriterios_impacto = tabla_subcriterios_impacto[!, ["Comparación"; subcriterios_impacto]]
end

# ╔═╡ c67e73d6-e475-4037-acf4-dccfa163cb5b
function viabilidad_pond_select(criterios::Dict)

	return combine() do Child
		
		inputs = [
			md""" $(v): $(
				Child(k, NumberField(0:0.05:100, default=0))
			)"""
			
			for (k, v) in criterios
		]
		
		md"""
		#### Comparación de subcriterios de Viabilidad: Indique cuan más importante es el primer subcriterio en comparación con el segundo.
		$(inputs)
		"""
	end
end

# ╔═╡ 7fd1b6fa-d125-4447-83c8-0682b4f432c3
@bind valor_viabilidad viabilidad_pond_select(Dict(
    :PP_vs_MO => "Prioridad Política versus Madurez Operativa"
))

# ╔═╡ 4a417cd0-c121-4998-9483-bedc724e3d55
begin 
	scv = valor_viabilidad.PP_vs_MO
	md""""""

end

# ╔═╡ 2ace2591-7f36-4447-8865-c6fa750ed45c
begin 
	matriz_subcriterios_viabilidad = [
    1   scv;
    1/scv 1
]
	md""""""

end

# ╔═╡ 4acecf52-9a9f-4abf-a914-9c335c25b521
begin
    tabla_subcriterios_viabilidad = DataFrame(matriz_subcriterios_viabilidad, Symbol.(subcriterios_viabilidad))

	# Metemos los nombres de los criterios como una nueva columna
    tabla_subcriterios_viabilidad[!, "Comparación"] = subcriterios_viabilidad

	# Nos aseguramos que "Comparación" sea la primera columna
    tabla_subcriterios_viabilidad = tabla_subcriterios_viabilidad[!, ["Comparación"; subcriterios_viabilidad]]
	
    tabla_subcriterios_viabilidad
end

# ╔═╡ 08420e09-65eb-4714-a6dc-967589481356
function soporte_pond_select(criterios::Dict)

	return combine() do Child
		
		inputs = [
			md""" $(v): $(
				Child(k, NumberField(0:0.05:100, default=0))
			)"""
			
			for (k, v) in criterios
		]
		
		md"""
		#### Comparación de subcriterios de Soporte Institucional: Indique cuan más importante es el primer subcriterio en comparación con el segundo.
		$(inputs)
		"""
	end
end

# ╔═╡ 4afb0b39-0fc6-464b-a65c-944fa501a50d
@bind valor_soporte soporte_pond_select(Dict(
    :VC_vs_MR => "Ventaja Comparativa del BID versus Movilización de Recursos"
))

# ╔═╡ c131fc4c-dfa7-49e3-a87c-c240f28610ff
begin 
	subcriterios_soporte_institucional = ["Ventaja Comparativa del BID", "Movilización de Recursos"]
	
	scs = valor_soporte.VC_vs_MR
	
	matriz_subcriterios_soporte_institucional = [
    1   scs;
    1/scs 1
]
	md""""""

end

# ╔═╡ 78007563-52d3-490f-915f-6c2271d82b4b
begin
    tabla_subcriterios_soporte_institucional = DataFrame(matriz_subcriterios_soporte_institucional, Symbol.(subcriterios_soporte_institucional))
	# Añadimos una columna con los nombres de los criterios
    tabla_subcriterios_soporte_institucional[!, "Comparación"] = subcriterios_soporte_institucional

    # Reordenamos las columnas para que "Comparación" sea la primera
    tabla_subcriterios_soporte_institucional = tabla_subcriterios_soporte_institucional[!, ["Comparación"; subcriterios_soporte_institucional]]
end

# ╔═╡ 182f6303-cecf-4601-a77a-784ff278c0dc
begin 
	md"""
	Matriz Principal: $(check_matrix(matriz_criterios_principales))
	
	Matriz Impacto: $(check_matrix(matriz_subcriterios_impacto))
	
	Matriz Viabilidad: $(check_matrix(matriz_subcriterios_viabilidad))
	
	Matriz Soporte Institucional: $(check_matrix(matriz_subcriterios_soporte_institucional))
	"""
end

# ╔═╡ 3f46e359-b88b-407a-9ca8-330da4b6e63d
begin 
	matrices = [
    (matriz_criterios_principales, "criterios principales"),
    (matriz_subcriterios_impacto, "subcriterios de impacto"),
    (matriz_subcriterios_viabilidad, "subcriterios de viabilidad"),
    (matriz_subcriterios_soporte_institucional, "subcriterios de soporte institucional"),
]	
	md""""""
end

# ╔═╡ 74ee7e98-4089-4c4c-b286-6afb48a4f7b1
for (matriz, nombre) in matrices
    CR = calcular_consistencia(matriz)
    if CR > 0.10
        println("Los valores de la matriz de $(nombre) deben ser corregidos.")
    else
        println("La matriz de $(nombre) es consistente.")
    end
end

# ╔═╡ 02b0a2ae-d908-4fa8-ae14-02b998937ac7
begin
	pesos_criterios_principales = calcular_pesos(matriz_criterios_principales)
	
	pesos_subcriterios_impacto = calcular_pesos(matriz_subcriterios_impacto)
	
	pesos_subcriterios_viabilidad = calcular_pesos(matriz_subcriterios_viabilidad)
	
	pesos_subcriterios_soporte_institucional = calcular_pesos(matriz_subcriterios_soporte_institucional)

	md""""""
	
end

# ╔═╡ c0fdba64-90fb-4440-99d9-9b72100c3168
begin
    # Crear vectores con los nombres de los criterios y los ponderadores
    nombres_criterios = criterios_principales
    ponderadores_criterios = pesos_criterios_principales

    nombres_subcriterios_impacto = subcriterios_impacto
    ponderadores_subcriterios_impacto = pesos_subcriterios_impacto

    nombres_subcriterios_viabilidad = subcriterios_viabilidad
    ponderadores_subcriterios_viabilidad = pesos_subcriterios_viabilidad

    nombres_subcriterios_soporte_institucional = subcriterios_soporte_institucional
    ponderadores_subcriterios_soporte_institucional = pesos_subcriterios_soporte_institucional

    # Creamos DataFrames con los resultados de los ponderadores y los nombres de los criterios
    df_criterios_principales = DataFrame(Criterio = nombres_criterios, Ponderador = ponderadores_criterios[:, 1])
    df_subcriterios_impacto = DataFrame(Subcriterio = nombres_subcriterios_impacto, Ponderador = ponderadores_subcriterios_impacto[:, 1])
    df_subcriterios_viabilidad = DataFrame(Subcriterio = nombres_subcriterios_viabilidad, Ponderador = ponderadores_subcriterios_viabilidad[:, 1])
    df_subcriterios_soporte_institucional = DataFrame(Subcriterio = nombres_subcriterios_soporte_institucional, Ponderador = ponderadores_subcriterios_soporte_institucional[:, 1])

    # Mostramos los DataFrames y vemos que anda todo bien
    df_criterios_principales, df_subcriterios_impacto, df_subcriterios_viabilidad, df_subcriterios_soporte_institucional
end

# ╔═╡ c7dfb9cb-3c72-4661-a3af-ae2c2c20979a
function criteria_selection(criteria::Vector)
	
	return combine() do Child
		
		inputs = [
			md""" $(name): $(
				Child(name, Slider(0:1:4))
			)"""
			
			for name in criteria
		]
		
		md"""
		#### Selecciona los valores para cada criterio
		$(inputs)
		"""
	end
end

# ╔═╡ 193417b9-6e7f-4696-9d71-203334c7e809
@bind project1 criteria_selection(["Impacto_Economico", "Impacto_Ambiental", "Impacto_Social", "Cooperacion", "Viabilidad_Prioridad_Política", "Viabilidad_Madurez_Operativa", "Soporte_Institucional_Ventaja_Comparativa_del_BID", "Soporte_Institucional_Movilizacion_de_Recursos"])

# ╔═╡ ee026067-6eef-4060-a5e6-39370e7338e5
valores_criterios = [
	project1.Impacto_Economico, 
	project1.Impacto_Ambiental, 
	project1.Impacto_Social, 
	project1.Cooperacion, 
	project1.Viabilidad_Prioridad_Política, 
	project1.Viabilidad_Madurez_Operativa, project1.Soporte_Institucional_Ventaja_Comparativa_del_BID, project1.Soporte_Institucional_Movilizacion_de_Recursos]

# ╔═╡ ad52b335-a8ef-4e0f-b9d3-3697adffba12
suma_impacto = sum(valores_criterios[1:3] .* pesos_subcriterios_impacto)

# ╔═╡ 27b30aa3-63ce-45a2-8409-10f358d063ce
suma_viabilidad = sum(valores_criterios[5:6] .* pesos_subcriterios_viabilidad)

# ╔═╡ a319f8cf-3f49-44a9-9bf8-d82476c2aebe
suma_soporte = sum(valores_criterios[7:8] .* pesos_subcriterios_soporte_institucional)

# ╔═╡ 389f081a-8e2d-4092-a75f-3d3f58b3ba07
valores_finales_criterios = [suma_impacto, valores_criterios[4], suma_viabilidad, suma_soporte]

# ╔═╡ e3875128-31ea-4415-8a1c-ac7eec8d632b
begin 
	valor_proyecto_1 = sum(valores_finales_criterios .* pesos_criterios_principales)

	md"""## El valor final del proyecto "Corredor pacífico de integración logística Mesoamericana" es: $(round(valor_proyecto_1, digits=2))"""

end

# ╔═╡ 9c047f19-9f28-4564-b465-726dd3a01248
@bind project2 criteria_selection(["Impacto_Economico", "Impacto_Ambiental", "Impacto_Social", "Cooperacion", "Viabilidad_Prioridad_Política", "Viabilidad_Madurez_Operativa", "Soporte_Institucional_Ventaja_Comparativa_del_BID", "Soporte_Institucional_Movilizacion_de_Recursos"])

# ╔═╡ 30504ae6-ec81-4762-93e6-0f10f16228ad
valores_criterios_2 = [
	project2.Impacto_Economico, 
	project2.Impacto_Ambiental, 
	project2.Impacto_Social, 
	project2.Cooperacion, 
	project2.Viabilidad_Prioridad_Política, 
	project2.Viabilidad_Madurez_Operativa, project2.Soporte_Institucional_Ventaja_Comparativa_del_BID, project2.Soporte_Institucional_Movilizacion_de_Recursos]

# ╔═╡ eb69a749-32ec-437b-97fa-897c4c2f185c
suma_impacto_2 = sum(valores_criterios_2[1:3] .* pesos_subcriterios_impacto)

# ╔═╡ 2bf051ae-0cb9-4694-afa6-698340145bb7
suma_viabilidad_2 = sum(valores_criterios_2[5:6] .* pesos_subcriterios_viabilidad)

# ╔═╡ a7b2cd32-d6dd-47df-b247-15d43b77531c
suma_soporte_2 = sum(valores_criterios_2[7:8] .* pesos_subcriterios_soporte_institucional)

# ╔═╡ b708e367-3a79-4abe-b1c1-775ccafee4a4
valores_finales_criterios_2 = [suma_impacto_2, valores_criterios_2[4], suma_viabilidad_2, suma_soporte_2]


# ╔═╡ 4f3baf8a-cabc-435b-8ffe-e0a3d1a8fbd8
begin 
	valor_proyecto_2 = sum(valores_finales_criterios_2 .* pesos_criterios_principales)
	
	md"""## El valor final del proyecto "Conglomerado portuario del Caribe" es: $(round(valor_proyecto_2, digits=2))"""

end

# ╔═╡ 5d5c24e1-e98d-44f5-9006-208cf3e2e078
@bind project3 criteria_selection(["Impacto_Economico", "Impacto_Ambiental", "Impacto_Social", "Cooperacion", "Viabilidad_Prioridad_Política", "Viabilidad_Madurez_Operativa", "Soporte_Institucional_Ventaja_Comparativa_del_BID", "Soporte_Institucional_Movilizacion_de_Recursos"])

# ╔═╡ 98f4ae84-f041-4938-9241-a84b485c47c0
valores_criterios_3 = [
	project3.Impacto_Economico, 
	project3.Impacto_Ambiental, 
	project3.Impacto_Social, 
	project3.Cooperacion, 
	project3.Viabilidad_Prioridad_Política, 
	project3.Viabilidad_Madurez_Operativa, project3.Soporte_Institucional_Ventaja_Comparativa_del_BID, project3.Soporte_Institucional_Movilizacion_de_Recursos]

# ╔═╡ 3989dfe8-5e39-4c2c-aa70-0b0f0c634e2e
suma_impacto_3 = sum(valores_criterios_3[1:3] .* pesos_subcriterios_impacto)

# ╔═╡ ebb752de-63a4-495b-8986-83c4f02f4e0d
suma_viabilidad_3 = sum(valores_criterios_3[5:6] .* pesos_subcriterios_viabilidad)

# ╔═╡ cf567fad-9ad4-4863-be8c-ecfe427d5f4e
suma_soporte_3 = sum(valores_criterios_3[7:8] .* pesos_subcriterios_soporte_institucional)

# ╔═╡ a2262bf2-d364-4960-b829-b12ac39d3489
valores_finales_criterios_3 = [suma_impacto_3, valores_criterios_3[4], suma_viabilidad_3, suma_soporte_3]

# ╔═╡ fcfe5013-6489-49d0-a71e-7bbdeeae716b
begin 
	valor_proyecto_3 = sum(valores_finales_criterios_3 .* pesos_criterios_principales)

	md"""## El valor final del proyecto "Integración eléctrica latinoamericana" es: $(round(valor_proyecto_3, digits=2))"""

end

# ╔═╡ 57442317-b142-4bce-ad1d-ef5626a009c5
@bind project4 criteria_selection(["Impacto_Economico", "Impacto_Ambiental", "Impacto_Social", "Cooperacion", "Viabilidad_Prioridad_Política", "Viabilidad_Madurez_Operativa", "Soporte_Institucional_Ventaja_Comparativa_del_BID", "Soporte_Institucional_Movilizacion_de_Recursos"])

# ╔═╡ 260f934d-97e9-4a7e-9a01-9f3fc45b3f94
valores_criterios_4 = [
	project4.Impacto_Economico, 
	project4.Impacto_Ambiental, 
	project4.Impacto_Social, 
	project4.Cooperacion, 
	project4.Viabilidad_Prioridad_Política, 
	project4.Viabilidad_Madurez_Operativa, project4.Soporte_Institucional_Ventaja_Comparativa_del_BID, project4.Soporte_Institucional_Movilizacion_de_Recursos]

# ╔═╡ a8a17f3f-a389-4ee4-95a4-0618612c8153
suma_impacto_4 = sum(valores_criterios_4[1:3] .* pesos_subcriterios_impacto)

# ╔═╡ 7733ddf1-fd04-41b1-844d-10b8f0992818
suma_viabilidad_4 = sum(valores_criterios_4[5:6] .* pesos_subcriterios_viabilidad)

# ╔═╡ 47e31c82-1f2a-4a8d-9ea6-6cbe1458652c
suma_soporte_4 = sum(valores_criterios_4[7:8] .* pesos_subcriterios_soporte_institucional)

# ╔═╡ 3f157169-f778-4e40-bda5-f7ab7eec193f
valores_finales_criterios_4 = [suma_impacto_4, valores_criterios_4[4], suma_viabilidad_4, suma_soporte_4]

# ╔═╡ 9d66192c-b62d-43dd-a97c-5d82de9bb2fe
begin 
	valor_proyecto_4 = sum(valores_finales_criterios_4 .* pesos_criterios_principales)

	md"""## El valor final del proyecto "Integración transfronteriza andina" es: $(round(valor_proyecto_4, digits=2))"""

end

# ╔═╡ 98018ba0-b92b-49a4-a0ed-abcade92d7cc
@bind project5 criteria_selection(["Impacto_Economico", "Impacto_Ambiental", "Impacto_Social", "Cooperacion", "Viabilidad_Prioridad_Política", "Viabilidad_Madurez_Operativa", "Soporte_Institucional_Ventaja_Comparativa_del_BID", "Soporte_Institucional_Movilizacion_de_Recursos"])

# ╔═╡ 5bdd7098-8c29-4269-ab21-51f7ed823539
valores_criterios_5 = [
	project5.Impacto_Economico, 
	project5.Impacto_Ambiental, 
	project5.Impacto_Social, 
	project5.Cooperacion, 
	project5.Viabilidad_Prioridad_Política, 
	project5.Viabilidad_Madurez_Operativa, project5.Soporte_Institucional_Ventaja_Comparativa_del_BID, project5.Soporte_Institucional_Movilizacion_de_Recursos]

# ╔═╡ 8cc58f8c-7957-40a2-9d83-353de6f76c11
suma_impacto_5 = sum(valores_criterios_5[1:3] .* pesos_subcriterios_impacto)

# ╔═╡ 023eac4d-d40e-43cb-bb3e-d5642d464599
suma_viabilidad_5 = sum(valores_criterios_5[5:6] .* pesos_subcriterios_viabilidad)

# ╔═╡ 00feb4c4-92ac-4b1d-be5f-9d1b5dbced3d
suma_soporte_5 = sum(valores_criterios_5[7:8] .* pesos_subcriterios_soporte_institucional)

# ╔═╡ b9c71c0a-706f-44f6-a9e2-9166ac4ab2c7
valores_finales_criterios_5 = [suma_impacto_5, valores_criterios_5[4], suma_viabilidad_5, suma_soporte_5]

# ╔═╡ 8deb8f4c-bfb9-47f6-bfd5-908259e00f5d
begin
	valor_proyecto_5 = sum(valores_finales_criterios_5 .* pesos_criterios_principales)

	md"""## El valor final del proyecto "Corredor bioceánico de capricornio" es: $(round(valor_proyecto_5, digits=2))"""

end

# ╔═╡ 8108f771-3a15-463c-bc5c-e01fa6600a94
@bind project6 criteria_selection(["Impacto_Economico", "Impacto_Ambiental", "Impacto_Social", "Cooperacion", "Viabilidad_Prioridad_Política", "Viabilidad_Madurez_Operativa", "Soporte_Institucional_Ventaja_Comparativa_del_BID", "Soporte_Institucional_Movilizacion_de_Recursos"])

# ╔═╡ ff15fbf3-df90-4f41-91dc-5ea219f20d33
valores_criterios_6 = [
	project6.Impacto_Economico, 
	project6.Impacto_Ambiental, 
	project6.Impacto_Social, 
	project6.Cooperacion, 
	project6.Viabilidad_Prioridad_Política, 
	project6.Viabilidad_Madurez_Operativa, project6.Soporte_Institucional_Ventaja_Comparativa_del_BID, project6.Soporte_Institucional_Movilizacion_de_Recursos]

# ╔═╡ fe0eb400-0b17-4b3a-9579-5b2c7420bbd4
suma_impacto_6 = sum(valores_criterios_6[1:3] .* pesos_subcriterios_impacto)

# ╔═╡ 01d66d3a-4947-43f3-9e17-050e69895d12
suma_viabilidad_6 = sum(valores_criterios_6[5:6] .* pesos_subcriterios_viabilidad)

# ╔═╡ 0186f246-7881-4d0f-a1b9-ef8be5b6e6e0
suma_soporte_6 = sum(valores_criterios_6[7:8] .* pesos_subcriterios_soporte_institucional)


# ╔═╡ b5d1b2fe-6c7b-4d13-9292-5042e17deb17
valores_finales_criterios_6 = [suma_impacto_6, valores_criterios_6[4], suma_viabilidad_6, suma_soporte_6]

# ╔═╡ f40b6cb6-6bac-483c-b643-dc22161653f7
begin
	valor_proyecto_6 = sum(valores_finales_criterios_6 .* pesos_criterios_principales)

	md"""## El valor final del proyecto "Iniciativa de la Cuenca del Plata" es: $(round(valor_proyecto_6, digits=2))"""
end


# ╔═╡ b62567f9-2338-4d01-b079-dc1355dd7b3b
@bind project7 criteria_selection(["Impacto_Economico", "Impacto_Ambiental", "Impacto_Social", "Cooperacion", "Viabilidad_Prioridad_Política", "Viabilidad_Madurez_Operativa", "Soporte_Institucional_Ventaja_Comparativa_del_BID", "Soporte_Institucional_Movilizacion_de_Recursos"])

# ╔═╡ 969cb1a9-ca20-45d1-b850-7aa2c8d5dfbd
valores_criterios_7 = [
	project7.Impacto_Economico, 
	project7.Impacto_Ambiental, 
	project7.Impacto_Social, 
	project7.Cooperacion, 
	project7.Viabilidad_Prioridad_Política, 
	project7.Viabilidad_Madurez_Operativa, project7.Soporte_Institucional_Ventaja_Comparativa_del_BID, project7.Soporte_Institucional_Movilizacion_de_Recursos
]

# ╔═╡ 8caedbb6-406b-4c74-b2fc-36143d3de488
suma_impacto_7 = sum(valores_criterios_7[1:3] .* pesos_subcriterios_impacto)

# ╔═╡ 9420fad8-914e-4e01-95b3-622dbfc6a2dd
suma_viabilidad_7 = sum(valores_criterios_7[5:6] .* pesos_subcriterios_viabilidad)

# ╔═╡ f269d306-47b2-4942-8f2e-b752e8725caf
suma_soporte_7 = sum(valores_criterios_7[7:8] .* pesos_subcriterios_soporte_institucional)

# ╔═╡ 8746fbdc-6960-4557-b8e9-0f699b61f24d
valores_finales_criterios_7 = [suma_impacto_7, valores_criterios_7[4], suma_viabilidad_7, suma_soporte_7]

# ╔═╡ a3633fc7-6d79-4135-a78e-88c8ec66bd6e
begin
	valor_proyecto_7 = sum(valores_finales_criterios_7 .* pesos_criterios_principales)

	md"""## El valor final del proyecto "Corredor de conectividad MERCOSUR-Chile" es: $(round(valor_proyecto_7, digits=2))"""
end

# ╔═╡ 73cded5d-0198-49e4-884b-f3a76e8bf397
@bind project8 criteria_selection(["Impacto_Economico", "Impacto_Ambiental", "Impacto_Social", "Cooperacion", "Viabilidad_Prioridad_Política", "Viabilidad_Madurez_Operativa", "Soporte_Institucional_Ventaja_Comparativa_del_BID", "Soporte_Institucional_Movilizacion_de_Recursos"])

# ╔═╡ d82e69e6-1662-4d92-93b0-ba923369765f
md"""Armamos un vector con todos los valores asignados a cada criterio"""

# ╔═╡ 33d40387-8404-402e-9980-7140f56c5abf
valores_criterios_8 = [
	project8.Impacto_Economico, 
	project8.Impacto_Ambiental, 
	project8.Impacto_Social, 
	project8.Cooperacion, 
	project8.Viabilidad_Prioridad_Política, 
	project8.Viabilidad_Madurez_Operativa, project8.Soporte_Institucional_Ventaja_Comparativa_del_BID, project8.Soporte_Institucional_Movilizacion_de_Recursos
]

# ╔═╡ fcf48ba5-ca20-4af8-8f3d-9f7c5fc29d58
md"""Calculamos la suma ponderada para las secciones con subcriterios"""

# ╔═╡ 18a7444c-4e20-4b41-98e7-021e337288cd
suma_impacto_8 = sum(valores_criterios_8[1:3] .* pesos_subcriterios_impacto)

# ╔═╡ 7e54e173-c053-46ac-af8e-9a87467f1f9a
suma_viabilidad_8 = sum(valores_criterios_8[5:6] .* pesos_subcriterios_viabilidad)

# ╔═╡ 80e056ed-763d-4cc1-872a-ed55238b68f0
suma_soporte_8 = sum(valores_criterios_8[7:8] .* pesos_subcriterios_soporte_institucional)

# ╔═╡ 727580c8-a680-4cce-b9b6-7b60e7249cde
valores_finales_criterios_8 = [suma_impacto_8, valores_criterios_8[4], suma_viabilidad_8, suma_soporte_8]

# ╔═╡ 92ccb385-b4fd-45de-a6d9-e2d65e7f8193
begin 
	valor_proyecto_8 = sum(valores_finales_criterios_8 .* pesos_criterios_principales)

	md"""## El valor final del proyecto "Conectividad bioceánica en la Patagonia" es: $(round(valor_proyecto_8, digits=2))"""

end

# ╔═╡ 3320a5c2-692c-465f-bb10-45a1c06fd099
md"""# Criterios para el proyecto "Infraestructura digital y de datos para la integración social y productiva en ALC"
"""

# ╔═╡ 4bcca1dc-f9f6-4e9e-86e5-8d0f5e8eb7d1
@bind project9 criteria_selection(["Impacto_Economico", "Impacto_Ambiental", "Impacto_Social", "Cooperacion", "Viabilidad_Prioridad_Política", "Viabilidad_Madurez_Operativa", "Soporte_Institucional_Ventaja_Comparativa_del_BID", "Soporte_Institucional_Movilizacion_de_Recursos"])

# ╔═╡ 02cffa89-a87f-4052-8190-ba44247b11e7
valores_criterios_9 = [
	project9.Impacto_Economico, 
	project9.Impacto_Ambiental, 
	project9.Impacto_Social, 
	project9.Cooperacion, 
	project9.Viabilidad_Prioridad_Política, 
	project9.Viabilidad_Madurez_Operativa, project9.Soporte_Institucional_Ventaja_Comparativa_del_BID, project9.Soporte_Institucional_Movilizacion_de_Recursos
]

# ╔═╡ d78f57cf-c00f-41c5-8a9a-d65f60d1d86f
suma_impacto_9 = sum(valores_criterios_9[1:3] .* pesos_subcriterios_impacto)

# ╔═╡ d06c93e8-590b-44e7-8037-49ae0e6c8c04
suma_viabilidad_9 = sum(valores_criterios_9[5:6] .* pesos_subcriterios_viabilidad)

# ╔═╡ 9df6f5ae-f33f-4bb9-a3de-9c7536bc6f7d
suma_soporte_9 = sum(valores_criterios_9[7:8] .* pesos_subcriterios_soporte_institucional)

# ╔═╡ a4ec57bf-0095-435a-9bc7-44a87cc856ec
valores_finales_criterios_9 = [suma_impacto_9, valores_criterios_9[4], suma_viabilidad_9, suma_soporte_9]

# ╔═╡ 455489af-ccbe-4779-863c-4301837d249e
begin 
	valor_proyecto_9 = sum(valores_finales_criterios_9 .* pesos_criterios_principales)

	md"""## El valor final del proyecto "Infraestructura digital y de datos para la integración social y productiva en ALC" es: $(round(valor_proyecto_9, digits=2))"""

end

# ╔═╡ f08fbff9-371c-472c-ace2-fed2b6465e62
md"""# Calculamos ahora el ranking final!"""

# ╔═╡ 8de5c2fe-408e-4063-9abf-155be1fec706
begin 

	valores_finales_proyectos = [valor_proyecto_1, valor_proyecto_2, valor_proyecto_3, valor_proyecto_4, valor_proyecto_5, valor_proyecto_6, valor_proyecto_7, valor_proyecto_8, valor_proyecto_9]

	nombres_proyectos = [
    "Corredor pacífico de integración logística Mesoamericana",
    "Conglomerado portuario del Caribe",
    "Integración eléctrica latinoamericana",
    "Integración transfronteriza andina",
    "Corredor bioceánico de capricornio",
    "Iniciativa de la Cuenca del Plata",
    "Corredor de conectividad MERCOSUR-Chile",
    "Conectividad bioceánica en la Patagonia",
	"Infraestructura digital y de datos para la integración social y productiva en ALC"]
	
	ranking_proyectos = sort([(valor, nombre) for (valor, nombre) in zip(valores_finales_proyectos, nombres_proyectos)], rev=true)
	
	md"""### Los valores finales para los proyectos son:

	$(nombres_proyectos[1]): $(round(valor_proyecto_1, digits = 2))

	$(nombres_proyectos[2]): $(round(valor_proyecto_2, digits = 2))

	$(nombres_proyectos[3]): $(round(valor_proyecto_3, digits = 2))

	$(nombres_proyectos[4]): $(round(valor_proyecto_4, digits = 2))

	$(nombres_proyectos[5]): $(round(valor_proyecto_5, digits = 2))

	$(nombres_proyectos[6]): $(round(valor_proyecto_6, digits = 2))

	$(nombres_proyectos[7]): $(round(valor_proyecto_7, digits = 2))

	$(nombres_proyectos[8]): $(round(valor_proyecto_8, digits = 2))

	$(nombres_proyectos[9]): $(round(valor_proyecto_9, digits = 2))
	"""
end

# ╔═╡ 39cad7f3-ac28-430a-9a90-cfbd92ab77ce
md"""# El ranking final queda conformado como:"""

# ╔═╡ a651aaba-2c1b-46b4-b27b-d1c36891f2b2
# Creando un DataFrame a partir de los rankings
ranking_df = DataFrame(
    Ranking = 1:length(ranking_proyectos),
    Nombre = [nombre for (valor, nombre) in ranking_proyectos],
    Valor = round.([valor for (valor, nombre) in ranking_proyectos], digits=2)
)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Cairo = "159f3aea-2a34-519c-b102-8c37f9878175"
Compose = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"
Images = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
Interact = "c601a237-2ae4-5e1e-952c-7a85b0c7eef1"
InteractiveUtils = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlotlyJS = "f0f68f2c-4968-5e81-91da-67840de0976a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
Pluto = "c3e4b0f8-55cb-11ea-2926-15256bba5781"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"
WebIO = "0f1e0344-ec1d-5b48-a673-e5cf874b6c29"

[compat]
Cairo = "~1.0.5"
Compose = "~0.9.5"
DataFrames = "~1.6.1"
FileIO = "~1.16.1"
HTTP = "~1.9.14"
Images = "~0.26.0"
Interact = "~0.10.5"
PlotlyJS = "~0.18.10"
Plots = "~1.38.17"
Pluto = "~0.19.27"
PlutoUI = "~0.7.52"
StatsPlots = "~0.15.6"
WebIO = "~0.8.21"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.1"
manifest_format = "2.0"
project_hash = "594aa5ec4e55860f81c8da5a3149b4807a31622c"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "91bd53c39b9cbfb5ef4b015e8b582d344532bd0a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.0"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "76289dc51920fdc6e0013c872ba9551d54961c24"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.2"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "9b9b347613394885fd1c8c7729bfc60528faa436"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.4"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "5ba6c757e8feccf03a1554dfaf3e26b3cfc7fd5e"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.1+1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "f83ec24f76d4c8f525099b2ac475fc098138ec31"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.4.11"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e5f08b5689b1aad068e01751889f2f615c7db36d"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.29"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AssetRegistry]]
deps = ["Distributed", "JSON", "Pidfile", "SHA", "Test"]
git-tree-sha1 = "b25e88db7944f98789130d7b503276bc34bc098e"
uuid = "bf4720bc-e11a-5d0c-854e-bdca1663c893"
version = "0.1.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "0c5f81f47bbbcf4aea7b2959135713459170798b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.5"

[[deps.Blink]]
deps = ["Base64", "Distributed", "HTTP", "JSExpr", "JSON", "Lazy", "Logging", "MacroTools", "Mustache", "Mux", "Pkg", "Reexport", "Sockets", "WebIO"]
git-tree-sha1 = "b1c61fd7e757c7e5ca6521ef41df8d929f41e3af"
uuid = "ad839575-38b3-5650-b840-f874b8c74a25"
version = "0.12.8"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "89e0654ed8c7aebad6d5ad235d6242c2d737a928"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.3"

[[deps.CSSUtil]]
deps = ["Colors", "JSON", "Markdown", "Measures", "WebIO"]
git-tree-sha1 = "b9fb4b464ec10e860abe251b91d4d049934f7399"
uuid = "70588ee8-6100-5070-97c1-3cb50ed05fe8"
version = "0.1.1"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays"]
git-tree-sha1 = "a0f80a09780eed9b1d106a1bf62041c2efc995bc"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e30f2f4e20f7f186dc36529910beaedc60cfa644"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.16.0"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "70232f82ffaab9dc52585e0dd043b5e0c6b714f1"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.12"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "b86ac2c5543660d238957dbde5ac04520ae977a7"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.15.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "02aa26a4cf76381be7f66e020a3eddeb27b0a092"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "dd3000d954d483c1aad05fe1eb9e6a715c97013e"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.22.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "e460f044ca8b99be31d35fe54fc33a5c33dd8ed7"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.9.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.Compose]]
deps = ["Base64", "Colors", "DataStructures", "Dates", "IterTools", "JSON", "LinearAlgebra", "Measures", "Printf", "Random", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "bf6570a34c850f99407b494757f5d7ad233a7257"
uuid = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
version = "0.9.5"

[[deps.ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "5372dbbf8f0bdb8c700db5367132925c0771ef7e"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.2.1"

[[deps.Configurations]]
deps = ["ExproniconLite", "OrderedCollections", "TOML"]
git-tree-sha1 = "434f446dbf89d08350e83bf57c0fc86f5d3ffd4e"
uuid = "5218b696-f38b-4ac9-8b61-a12ec717816d"
version = "0.17.5"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "f9d7112bfff8a19a3a4ea4e03a8e6a91fe8456bf"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.3"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.CustomUnitRanges]]
git-tree-sha1 = "1a3f97f907e6dd8983b744d2642651bb162a3f7a"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "1.0.2"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "04c738083f29f86e62c8afc341f0967d8717bdb8"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.6.1"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "b6def76ffad15143924a2199f72a5cd883a2e8a9"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.9"
weakdeps = ["SparseArrays"]

    [deps.Distances.extensions]
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "27a18994a5991b1d2e2af7833c4f8ecf9af6b9ea"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.99"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "e90caa41f5a86296e014e148ee061bd6c3edec96"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.ExproniconLite]]
deps = ["Pkg", "TOML"]
git-tree-sha1 = "d80b5d5990071086edf5de9018c6c69c83937004"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.3"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTViews]]
deps = ["CustomUnitRanges", "FFTW"]
git-tree-sha1 = "cbdf14d1e8c7c8aacbe8b19862e0179fd08321c2"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.3.2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "b4fbdd20c889804969571cc589900803edda16b7"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.7.1"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "299dc33549f68299137e51e6d49a13b5b1da9673"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "f372472e8672b1d993e93dada09e23139b509f9e"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.5.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.FunctionalCollections]]
deps = ["Test"]
git-tree-sha1 = "04cb9cfaa6ba5311973994fe3496ddec19b6292a"
uuid = "de31a74c-ac4f-5751-b3fd-e18cd04993ca"
version = "0.5.0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.FuzzyCompletions]]
deps = ["REPL"]
git-tree-sha1 = "e16dd964b4dfaebcded16b2af32f05e235b354be"
uuid = "fb4132e2-a121-4a70-b8a1-d5b831dcdcc2"
version = "0.5.1"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "8e2d86e06ceb4580110d9e716be26658effc5bfd"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.8"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "da121cbdc95b065da07fbb93638367737969693f"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.8+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "1cf1d7dcb4bc32d7b4a5add4232db3750c27ecb4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.8.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "cb56ccdd481c0dd7f975ad2b3b62d9eda088f7e2"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.9.14"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hiccup]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "6187bb2d5fcbb2007c39e7ac53308b0d371124bd"
uuid = "9fb69e20-1954-56bb-a84f-559cc56a8ff7"
version = "0.2.2"

[[deps.HistogramThresholding]]
deps = ["ImageBase", "LinearAlgebra", "MappedArrays"]
git-tree-sha1 = "7194dfbb2f8d945abdaf68fa9480a965d6661e69"
uuid = "2c695a8d-9458-5d45-9878-1b8a99cf7853"
version = "0.3.1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "d38bd0d9759e3c6cfa19bdccc314eccf8ce596cc"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.15"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "2e4520d67b0cef90865b3ef727594d2a58e0e1f8"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.11"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageBinarization]]
deps = ["HistogramThresholding", "ImageCore", "LinearAlgebra", "Polynomials", "Reexport", "Statistics"]
git-tree-sha1 = "f5356e7203c4a9954962e3757c08033f2efe578a"
uuid = "cbc4b850-ae4b-5111-9e64-df94c024a13d"
version = "0.3.0"

[[deps.ImageContrastAdjustment]]
deps = ["ImageBase", "ImageCore", "ImageTransformations", "Parameters"]
git-tree-sha1 = "eb3d4365a10e3f3ecb3b115e9d12db131d28a386"
uuid = "f332f351-ec65-5f6a-b3d1-319c6670881a"
version = "0.3.12"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "fc5d1d3443a124fde6e92d0260cd9e064eba69f8"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.1"

[[deps.ImageCorners]]
deps = ["ImageCore", "ImageFiltering", "PrecompileTools", "StaticArrays", "StatsBase"]
git-tree-sha1 = "24c52de051293745a9bad7d73497708954562b79"
uuid = "89d5987c-236e-4e32-acd0-25bd6bd87b70"
version = "0.1.3"

[[deps.ImageDistances]]
deps = ["Distances", "ImageCore", "ImageMorphology", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "08b0e6354b21ef5dd5e49026028e41831401aca8"
uuid = "51556ac3-7006-55f5-8cb3-34580c88182d"
version = "0.2.17"

[[deps.ImageFiltering]]
deps = ["CatIndices", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageBase", "ImageCore", "LinearAlgebra", "OffsetArrays", "PrecompileTools", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "c371a39622dc3b941ffd7c00e6b519d63b3f3f06"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.7.7"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "bca20b2f5d00c4fbc192c3212da8fa79f4688009"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.7"

[[deps.ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils"]
git-tree-sha1 = "b0b765ff0b4c3ee20ce6740d843be8dfce48487c"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.3.0"

[[deps.ImageMagick_jll]]
deps = ["JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "1c0a2295cca535fabaf2029062912591e9b61987"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "6.9.10-12+3"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "355e2b974f2e3212a75dfb60519de21361ad3cb7"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.9"

[[deps.ImageMorphology]]
deps = ["DataStructures", "ImageCore", "LinearAlgebra", "LoopVectorization", "OffsetArrays", "Requires", "TiledIteration"]
git-tree-sha1 = "6f0a801136cb9c229aebea0df296cdcd471dbcd1"
uuid = "787d08f9-d448-5407-9aad-5290dd7ab264"
version = "0.4.5"

[[deps.ImageQualityIndexes]]
deps = ["ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "LazyModules", "OffsetArrays", "PrecompileTools", "Statistics"]
git-tree-sha1 = "783b70725ed326340adf225be4889906c96b8fd1"
uuid = "2996bd0c-7a13-11e9-2da2-2f5ce47296a9"
version = "0.3.7"

[[deps.ImageSegmentation]]
deps = ["Clustering", "DataStructures", "Distances", "Graphs", "ImageCore", "ImageFiltering", "ImageMorphology", "LinearAlgebra", "MetaGraphs", "RegionTrees", "SimpleWeightedGraphs", "StaticArrays", "Statistics"]
git-tree-sha1 = "3ff0ca203501c3eedde3c6fa7fd76b703c336b5f"
uuid = "80713f31-8817-5129-9cf8-209ff8fb23e1"
version = "1.8.2"

[[deps.ImageShow]]
deps = ["Base64", "ColorSchemes", "FileIO", "ImageBase", "ImageCore", "OffsetArrays", "StackViews"]
git-tree-sha1 = "3b5344bcdbdc11ad58f3b1956709b5b9345355de"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.8"

[[deps.ImageTransformations]]
deps = ["AxisAlgorithms", "CoordinateTransformations", "ImageBase", "ImageCore", "Interpolations", "OffsetArrays", "Rotations", "StaticArrays"]
git-tree-sha1 = "7ec124670cbce8f9f0267ba703396960337e54b5"
uuid = "02fcd773-0e25-5acc-982a-7f6622650795"
version = "0.10.0"

[[deps.Images]]
deps = ["Base64", "FileIO", "Graphics", "ImageAxes", "ImageBase", "ImageBinarization", "ImageContrastAdjustment", "ImageCore", "ImageCorners", "ImageDistances", "ImageFiltering", "ImageIO", "ImageMagick", "ImageMetadata", "ImageMorphology", "ImageQualityIndexes", "ImageSegmentation", "ImageShow", "ImageTransformations", "IndirectArrays", "IntegralArrays", "Random", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "TiledIteration"]
git-tree-sha1 = "d438268ed7a665f8322572be0dabda83634d5f45"
uuid = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
version = "0.26.0"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3d09a9f60edf77f8a4d99f9e015e8fbf9989605d"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.7+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.IntegralArrays]]
deps = ["ColorTypes", "FixedPointNumbers", "IntervalSets"]
git-tree-sha1 = "be8e690c3973443bec584db3346ddc904d4884eb"
uuid = "1d092043-8f09-5a30-832f-7509e371ab51"
version = "0.1.5"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0cb9352ef2e01574eeebdb102948a58740dcaf83"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2023.1.0+0"

[[deps.Interact]]
deps = ["CSSUtil", "InteractBase", "JSON", "Knockout", "Observables", "OrderedCollections", "Reexport", "WebIO", "Widgets"]
git-tree-sha1 = "c5091992248c7134af7c90554305c600d5d9012b"
uuid = "c601a237-2ae4-5e1e-952c-7a85b0c7eef1"
version = "0.10.5"

[[deps.InteractBase]]
deps = ["Base64", "CSSUtil", "Colors", "Dates", "JSExpr", "JSON", "Knockout", "Observables", "OrderedCollections", "Random", "WebIO", "Widgets"]
git-tree-sha1 = "aa5daeff326db0a9126a225b58ca04ae12f57259"
uuid = "d3863d7c-f0c8-5437-a7b4-3ae773c01009"
version = "0.10.10"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.IntervalSets]]
deps = ["Dates", "Random"]
git-tree-sha1 = "8e59ea773deee525c99a8018409f64f19fb719e6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.7"
weakdeps = ["Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "4ced6667f9974fc5c5943fa5e2ef1ca43ea9e450"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.8.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "Requires", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "aa6ffef1fd85657f4999030c52eaeec22a279738"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.33"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSExpr]]
deps = ["JSON", "MacroTools", "Observables", "WebIO"]
git-tree-sha1 = "b413a73785b98474d8af24fd4c8a975e31df3658"
uuid = "97c1335a-c9c5-57fe-bc5d-ec35cebe8660"
version = "0.5.4"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "327713faef2a3e5c80f96bf38d1fa26f7a6ae29e"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.Kaleido_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "43032da5832754f58d14a91ffbe86d5f176acda9"
uuid = "f7e6163d-2fa5-5f23-b69c-1db539e41963"
version = "0.2.1+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "90442c50e202a5cdf21a7899c66b240fdef14035"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.7"

[[deps.Knockout]]
deps = ["JSExpr", "JSON", "Observables", "Test", "WebIO"]
git-tree-sha1 = "91835de56d816864f1c38fb5e3fad6eb1e741271"
uuid = "bcebb21b-c2e3-54f8-a781-646b90f6d2cc"
version = "0.2.6"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f689897ccbe049adb19a065c495e75f372ecd42b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.4+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "88b8f66b604da079a627b6fb2860d3704a6729a1"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.14"

[[deps.LazilyInitializedFields]]
git-tree-sha1 = "410fe4739a4b092f2ffe36fcb0dcc3ab12648ce1"
uuid = "0e77f7df-68c5-4e49-93ce-4cd80f5598bf"
version = "1.2.1"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "c3ce8e7420b3a6e071e0fe4745f5d4300e37b13f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.24"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "PrecompileTools", "SIMDTypes", "SLEEFPirates", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "c88a4afe1703d731b1c4fdf4e3c7e77e3b176ea2"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.165"

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.LoopVectorization.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "154d7aaa82d24db6d8f7e4ffcfe596f40bff214b"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2023.1.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.MetaGraphs]]
deps = ["Graphs", "JLD2", "Random"]
git-tree-sha1 = "1130dbe1d5276cb656f6e1094ce97466ed700e5a"
uuid = "626554b9-1ddb-594c-aa3c-2596fe9399a5"
version = "0.7.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.MsgPack]]
deps = ["Serialization"]
git-tree-sha1 = "fc8c15ca848b902015bd4a745d350f02cf791c2a"
uuid = "99f44e22-a591-53d1-9472-aa23ef4bd671"
version = "1.2.0"

[[deps.MultivariateStats]]
deps = ["Arpack", "LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "68bf5103e002c44adfd71fea6bd770b3f0586843"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.10.2"

[[deps.Mustache]]
deps = ["Printf", "Tables"]
git-tree-sha1 = "821e918c170ead5298ff84bffee41dd28929a681"
uuid = "ffc61752-8dc7-55ee-8c37-f3e9cdd09e70"
version = "1.0.17"

[[deps.Mux]]
deps = ["AssetRegistry", "Base64", "HTTP", "Hiccup", "MbedTLS", "Pkg", "Sockets"]
git-tree-sha1 = "0bdaa479939d2a1f85e2f93e38fbccfcb73175a5"
uuid = "a975b10e-0019-58db-a62f-e48ff68538c9"
version = "1.0.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "2c3726ceb3388917602169bed973dbc97f1b51a8"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.13"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "2ac17d29c523ce1cd38e27785a7d23024853a4bb"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.10"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "a4ca623df1ae99d09bc9868b008262d0c0ac1e4f"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "bbb5c2115d63c2f1451cb70e5ef75e8fe4707019"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.22+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "67eae2738d63117a196f497d7db789821bce61d1"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.17"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "9b02b27ac477cad98114584ff964e3052f656a0f"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.0"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "84a314e3926ba9ec66ac097e3635e270986b0f10"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.9+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.Pidfile]]
deps = ["FileWatching", "Test"]
git-tree-sha1 = "2d8aaf8ee10df53d0dfb9b8ee44ae7c04ced2b03"
uuid = "fa939f87-e72e-5be4-a000-7fc836dbe307"
version = "1.3.0"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f6cf8e7944e50901594838951729a1861e668cb8"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.2"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "56baf69781fc5e61607c3e46227ab17f7040ffa2"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.19"

[[deps.PlotlyJS]]
deps = ["Base64", "Blink", "DelimitedFiles", "JSExpr", "JSON", "Kaleido_jll", "Markdown", "Pkg", "PlotlyBase", "REPL", "Reexport", "Requires", "WebIO"]
git-tree-sha1 = "7452869933cd5af22f59557390674e8679ab2338"
uuid = "f0f68f2c-4968-5e81-91da-67840de0976a"
version = "0.18.10"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "9f8675a55b37a70aa23177ec110f6e3f4dd68466"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.17"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Pluto]]
deps = ["Base64", "Configurations", "Dates", "Distributed", "FileWatching", "FuzzyCompletions", "HTTP", "HypertextLiteral", "InteractiveUtils", "Logging", "LoggingExtras", "MIMEs", "Markdown", "MsgPack", "Pkg", "PrecompileSignatures", "PrecompileTools", "REPL", "RegistryInstances", "RelocatableFolders", "Sockets", "TOML", "Tables", "URIs", "UUIDs"]
git-tree-sha1 = "06fec2244568a4641e3352d20d0a0a608df6fa92"
uuid = "c3e4b0f8-55cb-11ea-2926-15256bba5781"
version = "0.19.27"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e47cd150dbe0443c3a3651bc5b9cbd5576ab75b7"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.52"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "240d7170f5ffdb285f9427b92333c3463bf65bf6"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.1"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase"]
git-tree-sha1 = "3aa2bb4982e575acd7583f01531f241af077b163"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "3.2.13"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.PrecompileSignatures]]
git-tree-sha1 = "18ef344185f25ee9d51d80e179f8dad33dc48eb1"
uuid = "91cefc8d-f054-46dc-8f8c-26e11d7c5411"
version = "3.0.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "9673d39decc5feece56ef3940e5dafba15ba0f81"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.1.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "ee094908d720185ddbdc58dbe0c1cbe35453ec7a"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "6ec7ac8412e83d57e313393220879ede1740f9ee"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.2"

[[deps.Quaternions]]
deps = ["LinearAlgebra", "Random", "RealDot"]
git-tree-sha1 = "da095158bdc8eaccb7890f9884048555ab771019"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.7.4"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegionTrees]]
deps = ["IterTools", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "4618ed0da7a251c7f92e869ae1a19c74a7d2a7f9"
uuid = "dee08c22-ab7f-5625-9660-a9af2021b33f"
version = "0.3.2"

[[deps.RegistryInstances]]
deps = ["LazilyInitializedFields", "Pkg", "TOML", "Tar"]
git-tree-sha1 = "ffd19052caf598b8653b99404058fce14828be51"
uuid = "2792f1a3-b283-48e8-9a74-f99dce5104f3"
version = "0.1.0"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays"]
git-tree-sha1 = "54ccb4dbab4b1f69beb255a2c0ca5f65a9c82f08"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.5.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "4b8586aece42bee682399c4c4aee95446aa5cd19"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.39"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "04bdff0b09c65ff3e06a05e3eb7b120223da3d39"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleWeightedGraphs]]
deps = ["Graphs", "LinearAlgebra", "Markdown", "SparseArrays"]
git-tree-sha1 = "4b33e0e081a825dbfaf314decf58fa47e53d6acb"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.4.0"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "c60ec5c62180f27efea3ba2908480f8055e17cee"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "7beb031cf8145577fbccacd94b8a8f4ce78428d3"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "f295e0a1da4ca425659c57441bcb59abb035a4bc"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.8"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "Requires", "SnoopPrecompile", "SparseArrays", "Static", "SuiteSparse"]
git-tree-sha1 = "33040351d2403b84afce74dae2e22d3f5b18edcb"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.4.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore"]
git-tree-sha1 = "9cabadf6e7cd2349b6cf49f1915ad2028d65e881"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.6.2"
weakdeps = ["Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "75ebe04c5bed70b91614d684259b661c9e6274a4"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.0"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StatsPlots]]
deps = ["AbstractFFTs", "Clustering", "DataStructures", "Distributions", "Interpolations", "KernelDensity", "LinearAlgebra", "MultivariateStats", "NaNMath", "Observables", "Plots", "RecipesBase", "RecipesPipeline", "Reexport", "StatsBase", "TableOperations", "Tables", "Widgets"]
git-tree-sha1 = "9115a29e6c2cf66cf213ccc17ffd61e27e743b24"
uuid = "f3b207a7-027a-5e70-b257-86293d7955fd"
version = "0.15.6"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "e383c87cf2a1dc41fa30c093b2a19877c83e1bc1"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.2.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "eda08f7e9818eb53661b3deb74e3159460dfbc27"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.2"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "8621f5c499a8aa4aa970b1ae381aae0ef1576966"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.4"

[[deps.TiledIteration]]
deps = ["OffsetArrays", "StaticArrayInterface"]
git-tree-sha1 = "1176cc31e867217b06928e2f140c90bd1bc88283"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.5.0"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "b7a5e99f24892b6824a954199a45e9ffcc1c70f0"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "64eb17acef1d9734cf09967539818f38093d9b35"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.16.2"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "b182207d4af54ac64cbc71797765068fdeff475d"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.64"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "ed8d92d9774b077c53e1da50fd81a36af3744c1c"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WebIO]]
deps = ["AssetRegistry", "Base64", "Distributed", "FunctionalCollections", "JSON", "Logging", "Observables", "Pkg", "Random", "Requires", "Sockets", "UUIDs", "WebSockets", "Widgets"]
git-tree-sha1 = "0eef0765186f7452e52236fa42ca8c9b3c11c6e3"
uuid = "0f1e0344-ec1d-5b48-a673-e5cf874b6c29"
version = "0.8.21"

[[deps.WebSockets]]
deps = ["Base64", "Dates", "HTTP", "Logging", "Sockets"]
git-tree-sha1 = "4162e95e05e79922e44b9952ccbc262832e4ad07"
uuid = "104b5d7c-a370-577a-8038-80a2059c5097"
version = "1.6.0"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "fcdae142c1cfc7d89de2d11e08721d0f2f86c98a"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.6"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╠═8d16ed10-40f1-11ee-22c7-0d2050858a0c
# ╟─e5ce8ca5-ac4c-437f-b6e7-0ad8df985358
# ╟─8aa54581-5d3b-42e4-aaf8-8d80d00b14d2
# ╟─f6cb3e82-5f54-4e3c-a13d-50eeac715463
# ╟─784194ab-5637-4458-9f6c-4d182c25ed41
# ╟─bf7490a5-400c-4885-9be3-323320e6c677
# ╟─d1a1cbbf-76ce-46bd-8a8f-f0dfb4ce996b
# ╟─29beb7d0-04f0-4620-b5c3-644583f066e4
# ╟─ff82c5c0-5139-4a4f-977c-f9193f639051
# ╟─81c4e439-7cf5-4874-9ce6-33f3fdd66515
# ╟─530978c0-5aae-46cd-b7c3-987d17442311
# ╟─281309e0-e252-4e06-bd49-c216e3ad0502
# ╟─ef636cd9-1d27-43fc-9323-cadc6a4df699
# ╟─795eecd5-8784-4dff-b753-36662af88dd3
# ╟─9de4a365-b0f1-46d4-991c-42188e4babcc
# ╟─eb74d652-7fe1-4cd3-9013-a11d3efb7f49
# ╟─7fd1b6fa-d125-4447-83c8-0682b4f432c3
# ╟─4a417cd0-c121-4998-9483-bedc724e3d55
# ╟─c67e73d6-e475-4037-acf4-dccfa163cb5b
# ╟─9a55981b-81e4-4601-89e5-0dbb41ce8f36
# ╟─2ace2591-7f36-4447-8865-c6fa750ed45c
# ╟─4acecf52-9a9f-4abf-a914-9c335c25b521
# ╟─4afb0b39-0fc6-464b-a65c-944fa501a50d
# ╟─08420e09-65eb-4714-a6dc-967589481356
# ╟─c131fc4c-dfa7-49e3-a87c-c240f28610ff
# ╟─78007563-52d3-490f-915f-6c2271d82b4b
# ╟─5a14d156-41be-47d8-84d7-e89aae7142a8
# ╟─9c116dd2-5e11-406c-83e4-1f3a7dd6d810
# ╟─182f6303-cecf-4601-a77a-784ff278c0dc
# ╟─1e8d91bd-ae61-4bcb-ab94-f262ebce1319
# ╟─6922f54f-4b48-4e0d-97c1-8469c4550f5c
# ╟─acec33c0-9012-42f2-91e7-f18dbb5de7e4
# ╟─3f46e359-b88b-407a-9ca8-330da4b6e63d
# ╟─74ee7e98-4089-4c4c-b286-6afb48a4f7b1
# ╟─8c007656-930c-4410-8d7c-7815d43b64c7
# ╟─6137955d-0201-4e57-a6e7-7dbb70b90645
# ╟─f716db17-1bdd-4d4e-9fe2-a578613bfb03
# ╟─ae3aa23c-0895-42ee-9989-7ef66459cdd1
# ╟─02b0a2ae-d908-4fa8-ae14-02b998937ac7
# ╟─c0fdba64-90fb-4440-99d9-9b72100c3168
# ╟─c938c87e-b87e-49db-a585-e3c6d80196f9
# ╟─193417b9-6e7f-4696-9d71-203334c7e809
# ╟─dfdf0b56-f9ad-427c-a600-5f56ce8596a0
# ╟─ee026067-6eef-4060-a5e6-39370e7338e5
# ╟─b0c81f63-8e86-4fa9-a764-5d9f7b772ea6
# ╟─ad52b335-a8ef-4e0f-b9d3-3697adffba12
# ╟─27b30aa3-63ce-45a2-8409-10f358d063ce
# ╟─a319f8cf-3f49-44a9-9bf8-d82476c2aebe
# ╟─2aab428f-f7ef-4a52-b928-13888e517db7
# ╟─389f081a-8e2d-4092-a75f-3d3f58b3ba07
# ╟─3b9f9560-fc6e-488a-a7fd-dcb71710f2e4
# ╟─e3875128-31ea-4415-8a1c-ac7eec8d632b
# ╟─a9e1f4a6-782b-44c2-9523-c9892731ccdf
# ╟─9c047f19-9f28-4564-b465-726dd3a01248
# ╟─c6e082bb-980d-4bb5-9e94-3d3197d6f699
# ╟─30504ae6-ec81-4762-93e6-0f10f16228ad
# ╟─62481b55-70f2-4525-811e-3d88a8bd41e1
# ╟─eb69a749-32ec-437b-97fa-897c4c2f185c
# ╟─2bf051ae-0cb9-4694-afa6-698340145bb7
# ╟─a7b2cd32-d6dd-47df-b247-15d43b77531c
# ╟─b708e367-3a79-4abe-b1c1-775ccafee4a4
# ╟─4f3baf8a-cabc-435b-8ffe-e0a3d1a8fbd8
# ╟─6043c869-030e-4728-bfe0-731498a8c9ad
# ╟─5d5c24e1-e98d-44f5-9006-208cf3e2e078
# ╟─98f4ae84-f041-4938-9241-a84b485c47c0
# ╟─08b56360-6a2f-4ee3-ab91-ad4d7c74a18b
# ╟─3989dfe8-5e39-4c2c-aa70-0b0f0c634e2e
# ╟─ebb752de-63a4-495b-8986-83c4f02f4e0d
# ╟─cf567fad-9ad4-4863-be8c-ecfe427d5f4e
# ╟─a2262bf2-d364-4960-b829-b12ac39d3489
# ╟─fcfe5013-6489-49d0-a71e-7bbdeeae716b
# ╟─7adc6966-4a50-414a-bbb5-38efc9ffd51c
# ╟─57442317-b142-4bce-ad1d-ef5626a009c5
# ╟─fa088abb-5d03-472e-b693-4eca3efc58b8
# ╟─260f934d-97e9-4a7e-9a01-9f3fc45b3f94
# ╟─93fa0b3e-3592-4140-aea3-237a414f9eaa
# ╟─a8a17f3f-a389-4ee4-95a4-0618612c8153
# ╟─7733ddf1-fd04-41b1-844d-10b8f0992818
# ╟─47e31c82-1f2a-4a8d-9ea6-6cbe1458652c
# ╟─3f157169-f778-4e40-bda5-f7ab7eec193f
# ╟─9d66192c-b62d-43dd-a97c-5d82de9bb2fe
# ╟─8c598b3c-f564-4d8f-b007-993f520b1589
# ╟─98018ba0-b92b-49a4-a0ed-abcade92d7cc
# ╟─5bdd7098-8c29-4269-ab21-51f7ed823539
# ╟─f75ffe35-159c-4e9e-8956-16e55816c094
# ╟─8cc58f8c-7957-40a2-9d83-353de6f76c11
# ╟─023eac4d-d40e-43cb-bb3e-d5642d464599
# ╟─00feb4c4-92ac-4b1d-be5f-9d1b5dbced3d
# ╟─b9c71c0a-706f-44f6-a9e2-9166ac4ab2c7
# ╟─8deb8f4c-bfb9-47f6-bfd5-908259e00f5d
# ╟─7e489c8d-e253-40ed-ba64-fb0672f5b3bb
# ╟─8108f771-3a15-463c-bc5c-e01fa6600a94
# ╟─ff15fbf3-df90-4f41-91dc-5ea219f20d33
# ╟─0e6ad1ac-a123-42dc-865a-d664ced846f0
# ╟─fe0eb400-0b17-4b3a-9579-5b2c7420bbd4
# ╟─01d66d3a-4947-43f3-9e17-050e69895d12
# ╟─0186f246-7881-4d0f-a1b9-ef8be5b6e6e0
# ╟─b5d1b2fe-6c7b-4d13-9292-5042e17deb17
# ╟─f40b6cb6-6bac-483c-b643-dc22161653f7
# ╟─9f2cc85b-b55e-4be5-8d68-b442a20b150d
# ╟─b62567f9-2338-4d01-b079-dc1355dd7b3b
# ╟─27aa5269-df7d-4d84-b671-32167ca254a8
# ╟─355aebc0-858d-4953-8a49-3e0315954e92
# ╟─969cb1a9-ca20-45d1-b850-7aa2c8d5dfbd
# ╟─ffd6ca37-0c54-49d7-8e0d-9f3f0601f0bf
# ╟─8caedbb6-406b-4c74-b2fc-36143d3de488
# ╟─9420fad8-914e-4e01-95b3-622dbfc6a2dd
# ╟─f269d306-47b2-4942-8f2e-b752e8725caf
# ╟─8746fbdc-6960-4557-b8e9-0f699b61f24d
# ╟─a3633fc7-6d79-4135-a78e-88c8ec66bd6e
# ╟─e863b703-5476-4297-87bb-aa378caf8070
# ╟─95768a04-cc66-422b-9a84-c67da26ed7cc
# ╟─73cded5d-0198-49e4-884b-f3a76e8bf397
# ╟─c7dfb9cb-3c72-4661-a3af-ae2c2c20979a
# ╟─d82e69e6-1662-4d92-93b0-ba923369765f
# ╟─33d40387-8404-402e-9980-7140f56c5abf
# ╟─fcf48ba5-ca20-4af8-8f3d-9f7c5fc29d58
# ╟─18a7444c-4e20-4b41-98e7-021e337288cd
# ╟─7e54e173-c053-46ac-af8e-9a87467f1f9a
# ╟─80e056ed-763d-4cc1-872a-ed55238b68f0
# ╟─727580c8-a680-4cce-b9b6-7b60e7249cde
# ╟─92ccb385-b4fd-45de-a6d9-e2d65e7f8193
# ╟─3320a5c2-692c-465f-bb10-45a1c06fd099
# ╟─4bcca1dc-f9f6-4e9e-86e5-8d0f5e8eb7d1
# ╟─02cffa89-a87f-4052-8190-ba44247b11e7
# ╟─d78f57cf-c00f-41c5-8a9a-d65f60d1d86f
# ╟─d06c93e8-590b-44e7-8037-49ae0e6c8c04
# ╟─9df6f5ae-f33f-4bb9-a3de-9c7536bc6f7d
# ╟─a4ec57bf-0095-435a-9bc7-44a87cc856ec
# ╟─455489af-ccbe-4779-863c-4301837d249e
# ╟─f08fbff9-371c-472c-ace2-fed2b6465e62
# ╟─8de5c2fe-408e-4063-9abf-155be1fec706
# ╟─39cad7f3-ac28-430a-9a90-cfbd92ab77ce
# ╟─a651aaba-2c1b-46b4-b27b-d1c36891f2b2
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
