from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt

#Notacion IUPAC de los nucleotidos (incluyendo las bases canónicas y las degeneradas)

diccionario_bases = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'}, 'U': {'T'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'C', 'G'},
    'W': {'A', 'T'}, 'K': {'G', 'T'}, 'M': {'A', 'C'},
    'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'T'},
}

#Diccionario que asigna el complementario a una base

diccionario_basesComplementarias = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
        'D': 'H', 'H': 'D', 'N': 'N',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
        'r': 'y', 'y': 'r', 's': 's', 'w': 'w',
        'k': 'm', 'm': 'k', 'b': 'v', 'v': 'b',
        'd': 'h', 'h': 'd', 'n': 'n'
    }


def compatibilidad(base1, base2):
    '''La entrada son dos bases, se van a buscar en el mismo diccionario (diccionario_bases) usando .get
    El valor asignado para base1, se guarda como un diccionario en la variable bases_set1
    El valor asignado para base2 se guarda como un diccionario en la variable bases_set2
    Se comparan los elementos de ambos diccionarios y si existe por lo menos uno en ambos (interseccion) retorna True
    Si la comparación arroja una lista vacía retorna False'''
    bases_set1 = diccionario_bases.get(base1, {base1})
    bases_set2 = diccionario_bases.get(base2, {base2})
    resultado_compatibilidad = bool(bases_set1.intersection(bases_set2))
    return resultado_compatibilidad


def limpiar_degenerados(secuencia):
    '''Convierte todas las bases en una secuencia en mayusculas
    Recorre la secuencia, base por base y si coindice con algun elemento del conjunto validos,
    lo agrega a una lista vacía llamada secuencia_limpia, si no, lo reemplaza por 'N'
    finalmente, retorna la secuencia_limpia convertida en un string'''
    validos = set('ACGTN')
    secuencia_limpia = []
    secuencia = secuencia.upper()
    for base in secuencia:
        if base in validos:
            secuencia_limpia.append(base)
        else:
            secuencia_limpia.append('N')
    secuencia_final = "".join(secuencia_limpia)
    return secuencia_final


def reverso_complementario(secuencia):
    '''Invierte la secuencia de entrada (str) y la agrega en una variable llamada secuencia_invertida
    Reemplaza las bases de la secuencia_invertida por su correspondiente base complementaria
    y traduce esa lista en un string'''
    secuencia_invertida = secuencia[::-1]
    lista_basesComplementarias = []
    for base in secuencia_invertida:
        baseComplementaria = diccionario_basesComplementarias.get(base, base)
        lista_basesComplementarias.append(baseComplementaria)
    secuencia_complementaria = "".join(lista_basesComplementarias)
    return secuencia_complementaria


def distancia_levenshtein(s1, s2):
    '''Calcula la Distancia de Levenshtein entre dos cadenas.
    La distancia es el número mínimo de operaciones (inserción, eliminación o sustitución) de un solo carácter
    necesarias para transformar la cadena s1 en la cadena s2.
    Esta implementación asume que el costo de Inserción, Eliminación y Sustitución es 1,
    excepto cuando los caracteres coinciden (costo 0).
    Usa la función compatibilidad'''
    numero_filas = len(s1) + 1
    numero_columnas = len(s2) + 1
    matriz_distancias = [[0] * numero_columnas for _ in range(numero_filas)] #crea la matriz y la llena con ceros

    for i in range(numero_filas):
        matriz_distancias[i][0] = i #Costo de eliminar i caracteres de s1 para llegar a una cadena vacía

    for j in range(numero_columnas):
        matriz_distancias[0][j] = j #Costo de insertar j caracteres en s1 (vacía) para llegar a s2[:j]

    for i in range(1, numero_filas):
        for j in range(1, numero_columnas):
            if compatibilidad(s1[i - 1], s2[j - 1]):
                costo_sustitucion = 0
            else:
                costo_sustitucion = 1

            costo_eliminacion = matriz_distancias[i - 1][j] + 1
            costo_insercion = matriz_distancias[i][j - 1] + 1
            costo_sustitucion = matriz_distancias[i - 1][j - 1] + costo_sustitucion
            matriz_distancias[i][j] = min(costo_eliminacion, costo_insercion, costo_sustitucion)
    return matriz_distancias[numero_filas - 1][numero_columnas - 1]


def smith_waterman_score(s1, s2, match=1, mismatch=-1, gap=-2):
    '''Calcula la máxima puntuación de alineamiento local entre dos secuencias (s1 y s2)
    Args:
        s1 (str): La primera secuencia (e.g., ADN o proteína).
        s2 (str): La segunda secuencia.
        match (int): Puntuación por una coincidencia (default: 1).
        mismatch (int): Puntuación por una no coincidencia (default: -1).
        gap (int): Penalización por abrir o extender un hueco (gap) (default: -2).

    Returns:
        int: La máxima puntuación de alineamiento local encontrada en la matriz.'''
    numero_filas = len(s1) + 1
    numero_columnas = len(s2) + 1
    matriz_puntuacion = [[0] * numero_columnas for _ in range(numero_filas)]
    max_score = 0 #maxima puntuacion encontrada, para el alineamiento local
    for i in range(1, numero_filas):
        for j in range(1, numero_columnas):

            #Caracteres que se están comparando
            caracter1 = s1[i - 1]
            caracter2 = s2[j - 1]

            #Movimiento diagonal (coincidencia o no coincidencia)
            puntuacion_sustitucion = match if caracter1 == caracter2 else mismatch
            puntuacion_diagonal = matriz_puntuacion[i - 1][j - 1] + puntuacion_sustitucion

            #Movimiento superior (Hueco en s2 / Eliminación en s1)
            puntuacion_arriba = matriz_puntuacion[i - 1][j] + gap

            #Movimiento izquierdo (Hueco en s1 / Inserción en s2)
            puntuacion_izquierda = matriz_puntuacion[i][j - 1] + gap

            puntuacion_actual = max(0, puntuacion_diagonal, puntuacion_arriba, puntuacion_izquierda)
            matriz_puntuacion[i][j] = puntuacion_actual
            if puntuacion_actual > max_score:
                max_score = puntuacion_actual
    return max_score


def buscar_mejor_match_levenshtein(primer, secuencia_objetivo):
    '''Busca la subcadena en la secuencia objetivo que tiene la menor Distancia de Levenshtein
    con respecto al primer de entrada.

    Simula una 'búsqueda local por ventana deslizante' para encontrar el sitio de unión más probable para un primer.

    Returns:
        tuple[int, float, str]: que contiene:
            1. mejor_posicion (int): La posición inicial del mejor match en la secuencia objetivo.
            2. mejor_distancia (float): La distancia de Levenshtein más baja encontrada.
            3. mejor_subcadena (str): La subcadena de la secuencia objetivo que produjo la mejor distancia.

    Retorna (-1, float('inf'), '') si el cebador es más largo que la secuencia objetivo.'''
    len_primer = len(primer)
    mejor_distancia = float('inf')
    mejor_subcadena = ""
    mejor_posicion = -1

    if len_primer > len(secuencia_objetivo):
        return -1, float('inf'), ""

    ancho_ventana = len_primer
    rango_busqueda = len(secuencia_objetivo) - len_primer + 1
    for i in range(rango_busqueda):
        subcadena_actual = secuencia_objetivo[i: i + ancho_ventana]
        distancia_actual = distancia_levenshtein(primer, subcadena_actual)
        if distancia_actual < mejor_distancia:
            mejor_distancia = distancia_actual
            mejor_subcadena = subcadena_actual
            mejor_posicion = i
            if distancia_actual == 0:
                break
    return mejor_posicion, mejor_distancia, mejor_subcadena


def generar_string_alineamiento(s1, s2):
    '''Genera una cadena de caracteres que representa visualmente el alineamiento entre dos secuencias de igual longitud.

    Args:
        secuencia_a (str): La primera secuencia alineada.
        secuencia_b (str): La segunda secuencia alineada.

    Returns:
        str: Una cadena de la misma longitud que las entradas, donde:
            - '|' (Barra vertical) indica una coincidencia (match).
            - ' ' (Espacio) indica una no coincidencia o sustitución (mismatch)'''
    alineamiento = []
    for i in range(len(s1)):
        caracter1 = s1[i]
        caracter2 = s2[i]

        if compatibilidad(caracter1, caracter2):
            alineamiento.append('|')
        else:
            alineamiento.append(' ')
    return "".join(alineamiento)


def parametros_primer(primer_fwd, primer_rev):
    '''Calcula GC, Temperaturas de Fusión (Tm) y autocomplementariedad (tendencia a formar dímeros) para un par de primers
        Usa fórmulas empíricas para la Tm y la distancia de Levenshtein para la autocomplementariedad.

    Args:
        cebador_fwd (str): Secuencia del cebador Forward (hacia adelante).
        cebador_rev (str): Secuencia del cebador Reverse (hacia atrás).

    Returns:
        tuple[float, float, float, float, int, int]: Una tupla con los resultados redondeados:
            1. porc_GC_fwd (float): Porcentaje de GC del cebador Forward.
            2. porc_GC_rev (float): Porcentaje de GC del cebador Reverse.
            3. Tm_fwd (float): Temperatura de Fusión del cebador Forward (°C).
            4. Tm_rev (float): Temperatura de Fusión del cebador Reverse (°C).
            5. autocomplementariedad_fwd (int): Distancia de Levenshtein entre el cebador fwd y su reverso complementario.
            6. autocomplementariedad_rev (int): Distancia de Levenshtein entre el cebador rev y su reverso complementario'''
    ##Calculo de GC
    #Primer forward (fwd)
    G_fwd_count = primer_fwd.count('G')
    C_fwd_count = primer_fwd.count('C')
    GC_fwd = G_fwd_count + C_fwd_count
    GC_fwd_porc = (GC_fwd / len(primer_fwd)) * 100

    #Primer reverse (rev)
    G_rev_count = primer_rev.count('G')
    C_rev_count = primer_rev.count('C')
    GC_rev = G_rev_count + C_rev_count
    GC_rev_porc = (GC_rev / len(primer_rev)) * 100

    ##Calculo de Tm
    Tm_fwd = 64.9 + (41 * (GC_fwd - 16.4)) / len(primer_fwd)
    Tm_rev = 64.9 + (41 * (GC_rev - 16.4)) / len(primer_rev)

    ##Calculo de autocomplementariedad
    #Primer forward (fwd)
    rc_fwd = reverso_complementario(primer_fwd)
    autocomplemetary_fwd = distancia_levenshtein(primer_fwd, rc_fwd)

    #Primer reverse (rev)
    rc_rev = reverso_complementario(primer_rev)
    autocomplemetary_rev = distancia_levenshtein(primer_rev, rc_rev)

    return (round(GC_fwd_porc, 2),
            round(GC_rev_porc, 2),
            round(Tm_fwd, 2),
            round(Tm_rev, 2),
            autocomplemetary_fwd,
            autocomplemetary_rev)

def crear_mapa_genetico(longitud_total, lista_cds, lista_amplicones):
    '''Genera la figura de Matplotlib del mapa genético que contiene la ubicacion de los CDS y los amplicones
    Usa la libreria DnaFeaturesViewer
    Args:
        longitud_total (int): La longitud total de la secuencia a mapear.
        lista_cds (list): Una lista de diccionarios, donde cada diccionario describe un CDS
        con claves como "start", "end", "strand" y "name".
        lista_amplicones (list): Una lista de diccionarios que describen regiones de amplificación
        con claves como "start", "end" y "name".

    Returns:
        el objeto 'fig' para Streamlit'''
    caracteristicas_graficas = []

    # Procesar y agregar los CDS (AZUL)
    for cds in lista_cds:
        if str(cds['strand']) in ['1', '+', 'Derecha']:
            direccion_hebra = 1
        else:
            direccion_hebra = -1

        caracteristicas_graficas.append(GraphicFeature(start=cds["start"], end=cds["end"], strand=direccion_hebra, color="#2c7fb8", label=cds["name"]))

    # Procesar y agregar los Amplicones (ROJO/NARANJA)
    for amplicon in lista_amplicones:
        caracteristicas_graficas.append(GraphicFeature(start=amplicon["start"], end=amplicon["end"], strand=None, color="#e34a33", label=amplicon["name"]))

    # Crear el objeto de registro gráfico
    registro_grafico = GraphicRecord(sequence_length=longitud_total, features=caracteristicas_graficas)

    # Dibujar y configurar la figura (fig)
    figura, eje_plot = plt.subplots(1, 1, figsize=(15, 3))
    registro_grafico.plot(ax=eje_plot, figure_width=15)
    eje_plot.set_xlim(0, longitud_total)

    return figura

def preparar_entradas(string_objetivo):
    string_limpio = string_objetivo.strip().upper()
    return string_limpio


def obtener_veredictoFinal(fwd_upper: str, subcadena_fwd: str,
                            rc_rev_upper: str, subcadena_rev: str) -> bool:
    """
    Evalúa si un par de cebadores (Forward y Reverse) cumple con los criterios
    de calidad de alineamiento para la PCR, basándose en la Distancia de Levenshtein.

    Criterios de validación:
        1. Distancia Total (Match global): <= 5 ediciones permitidas.
        2. Distancia Anclaje 3' (Match perfecto en el extremo de iniciación): 0 ediciones requeridas.

    Args:
        fwd_upper (str): Secuencia del cebador Forward.
        subcadena_fwd (str): Subcadena objetivo alineada con el cebador FWD.
        rc_rev_upper (str): Secuencia del cebador Reverse Complementario.
        subcadena_rev (str): Subcadena objetivo alineada con el cebador REV.

    Returns:
        bool: True si AMBOS cebadores pasan TODOS los criterios de validación; False en caso contrario.
    """

    # --- 1. CÁLCULO DE LA DISTANCIA TOTAL (Match Global) ---

    # Se calcula la distancia de Levenshtein para el cebador completo.
    dist_fwd = distancia_levenshtein(fwd_upper, subcadena_fwd)
    dist_rev = distancia_levenshtein(rc_rev_upper, subcadena_rev)

    # --- 2. CÁLCULO DE LA DISTANCIA DE ANCLAJE 3' (Extremo Crítico) ---

    # Se calcula la distancia de Levenshtein para los últimos 5 nucleótidos (extremo 3').
    # Un match perfecto (distancia 0) es crucial aquí.

    # FWD: Comparar los últimos 5 caracteres del cebador y la subcadena
    dist_3_fwd = distancia_levenshtein(fwd_upper[-5:], subcadena_fwd[-5:])

    # REV: Comparar los últimos 5 caracteres del cebador y la subcadena
    dist_3_rev = distancia_levenshtein(rc_rev_upper[-5:], subcadena_rev[-5:])

    # --- 3. EVALUACIÓN DE CRITERIOS (Pasa/Falla) ---

    # Criterio 1: La distancia total debe ser tolerablemente baja (<= 5).
    pasa_dist_total_fwd = (dist_fwd <= 5)
    pasa_dist_total_rev = (dist_rev <= 5)

    # Criterio 2: El anclaje 3' debe ser un match perfecto (== 0).
    pasa_anclaje_3_fwd = (dist_3_fwd == 0)
    pasa_anclaje_3_rev = (dist_3_rev == 0)

    # --- 4. VEREDICTO FINAL ---

    # El par de cebadores pasa si todas las condiciones son True.
    veredicto_final = (pasa_dist_total_fwd and pasa_dist_total_rev and
                       pasa_anclaje_3_fwd and pasa_anclaje_3_rev)

    return veredicto_final