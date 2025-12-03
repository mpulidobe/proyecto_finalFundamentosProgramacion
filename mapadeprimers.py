import numpy as np
import streamlit as st
import primer3
from electroforesis_gel import simular_electroforesis
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord
import io

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


# STREAMLIT

st.set_page_config(page_title="GenVis", layout="wide")
st.title("GenVis - Herramienta Bioinformática para el Análisis y Visualización de Secuencias")

# ENTRADA DE SECUENCIA PRINCIPAL (Siempre visible)
st.header("Entrada de la Secuencia Objetivo")
default_sequence = """GTCATCTTTTATTCTTAATCAAACCTCACTCAGAAAAATCCAGAACTGTNNNAATAAGAACCAAAGCCACAVSRATGTCTACAAGTGTTTATGAATCGATCATTCAGACGAAAGCTTCGGTCTGGGGATCTACTGCATCTGGCAAATCTATTGTGGACTCTTACTGGATCCACGAGTTTTTAACTGGTTCTCCATTGATTCAAACTCAGTTGTATTCTGATTCAAGAAGCAAAAGCAGCTTTGGCTACACCACACGAGTTGGTGATCTTCVDSBCTTCAGAAGAGAAAGAGATTCTTTCTCAACACTTGTACATCCCTATTTTTGATGACATTGATTTCAACATCAATATCAATGATTCAGTCATGACAGTATCCGTTTGCTCCAACACGGTCAATGCTAATGGAGTGAAACATCAGGGTCATCTGAAGGTGCTTTCTCTTGCTCAACTGCACTCTATAGAGCCTACAATGAGCCGATCTGACATTGCTGACAGATTCCGTCTTCAAGAAAAAGACGTGATTCCCAATGACAGATACATTGATGCTGCTAACAAAGGCTCTCTTTCATGTGTTAAAGAGCATTCCTATAAAGTCGAAATGTGCCACAACCAAGCATTAGGTAAAGTTAATGTTCTATCCCCTAACAGAAATGTTCATGAATGGCTGTACAGCTTCAAGCCAGCTTTCAACCAGATTGAAAGCAACAACAGAACTGTAAATTCTCTTGCAGTGAAATCTCTGCTCATGTCTGCAGAAAACAACATAATGCCTAACTCTCAGGCCTTTGTCAAAGCTTCTACAGGCTCCCAGTTCAAGCTAAACCTCTGGCTGAGGATTCCTAAAGTTCTGAAACAGGTTTCTATTCAGAAACTATTTAAAGTTGCAGAAGATGAAACAGACAAAAGCTTTTATTTGTCTATTGCTTGCATCCCTAACCACAACAGCGTTGAAACAGCVSTMCTTGAATGTGACCATCATCTGCAAGCATCAGCTCCCAATCTCGAAGCTCAAAGCCCCTTTWKDTGAATTGACAATGATGTTTTCTGATCTGAGAGAACCTTACAACGTTGTGCATGATCCTTCTTACCCTCAAAGAATTGTTCATGCTCTGCTTGAGACACACACATCTTTTGCCCAAACTCTTTGCAATAACTTGCAAGAAGATGTGGTCATCTACACTTTGAACAACCCTGAGCTGACTTCTTTAAAGTTAGATTTAGGTAAGAAAACCCTAAATTACAGTGAAGATGCTTATAATAAGAAATATTTTCTTTCAAAAACTCTTGAATGCCTCCCAGTAAACACACAGACTATGTCTTATTTAGACAGCATTCAAATTCCCTCATGGAAGATTGACTTTGCCAGAGGAGAAATCAAAATTTCCCCTCAATCAATCTCTGTTGCAAAATCTTTgttgaagctagatcttgatgtgatcagaggaaagaaatctctgcctcagggagcttctgaatcagagtcaaagcaatttgtgtctatttgtctgctcctttaactattcctttcctctttaaatcttctttcagcttctttccaactfctttcaatctcttttcargtttctctcatttcctttaaatttcctaattttctttacttcctttct"""

seq_input = st.text_area("A continuación ingrese la secuencia",
                         value=default_sequence,
                         height=250, key="seq_input_main")

st.divider()

# NAVEGACIÓN POR PESTAÑAS
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "🧬 Analizar Par de Primers (PCR)", # PESTAÑA 1
    "🔬 Generar Primers Nuevos", # PESTAÑA 2
    "📊 Analizar Múltiples Pares", # PESTAÑA 3
    "🗺️ Visualizador Mapa Genético",  # PESTAÑA 4
    "💻 Simular Electroforesis"  # PESTAÑA 5
])

# =======================================================================
# PESTAÑA 1: ANALIZAR PAR DE PRIMERS
# =======================================================================
with tab1:
    st.header("Analizar un Par de Primers (Forward y Reverse)")
    st.info(
        "Esta pestaña aplica un filtro biológico: la distancia total debe ser <= 5 Y los últimos 5pb del extremo 3' deben tener una distancia == 0.")

    default_rev = "GGGGCTTTGAGCTTCGAGAT"
    default_fwd = "ACAACAGCGTTGAAACAGCC"

    col1, col2 = st.columns(2)
    with col1:
        fwd_input = st.text_input("Primer Forward", value=default_fwd, key="fwd_input_tab1")
    with col2:
        rev_input = st.text_input("Primer Reverse", value=default_rev, key="rev_input_tab1")

    if st.button("Analizar Par", use_container_width=True, type="primary", key="btn_analizar_par"):
        st.header("Resultados del Análisis del Par")

        if not seq_input or not fwd_input or not rev_input:
            st.warning("Por favor, ingresa la secuencia target y ambos primers")
        else:
            seq_upper = seq_input.strip().upper()
            fwd_upper = fwd_input.strip().upper()
            rev_upper = rev_input.strip().upper()

            if len(seq_upper) < 70:
                st.error(f"Error: El genoma de ({len(seq_upper)} pb) es más corto que el mínimo de 70 pb")
            else:
                st.subheader("Resultados de la Búsqueda")
                with st.spinner("Buscando el mejor match para ambos primers..."):
                    rc_rev_upper = reverso_complementario(rev_upper)
                    posicion_fwd, dist_fwd, subcadena_fwd = buscar_mejor_match_levenshtein(fwd_upper, seq_upper)
                    posicion_rev, dist_rev, subcadena_rev = buscar_mejor_match_levenshtein(rc_rev_upper, seq_upper)

                if posicion_fwd == -1:
                    st.error("Error: El primer Forward no se encontró en la secuencia")
                elif posicion_rev == -1:
                    st.error("Error: El Reverso Complementario del primer Reverse no se encontró en la secuencia")
                elif posicion_fwd >= posicion_rev:
                    st.error(
                        f"Error: El primer Forward (posición {posicion_fwd}) se encontró DESPUÉS del primer Reverse (posición {posicion_rev}).")
                else:
                    st.success("Primeras encontrados. Aplicando filtros biológicos...")

                    dist_3_fwd = distancia_levenshtein(fwd_upper[-5:], subcadena_fwd[-5:])
                    dist_3_rev = distancia_levenshtein(rc_rev_upper[-5:], subcadena_rev[-5:])

                    pasa_dist_total_fwd = (dist_fwd <= 5)
                    pasa_dist_total_rev = (dist_rev <= 5)
                    pasa_anclaje_3_fwd = (dist_3_fwd == 0)
                    pasa_anclaje_3_rev = (dist_3_rev == 0)

                    veredicto_final = (pasa_dist_total_fwd and pasa_dist_total_rev and
                                       pasa_anclaje_3_fwd and pasa_anclaje_3_rev)

                    if veredicto_final:
                        st.success(
                            "VEREDICTO: ¡ÉXITO! El par de primers CUMPLE con todas las reglas de anclaje y distancia.")
                    else:
                        st.error(
                            "VEREDICTO: FALLO. El par de primers NO cumple con los criterios de anclaje y/o distancia total. Revisa los detalles abajo.")

                    inicio_fwd = posicion_fwd
                    fin_rev = posicion_rev + len(rc_rev_upper)
                    tamano_amplicon = fin_rev - inicio_fwd
                    st.metric("Tamaño Estimado del Amplicón", f"{tamano_amplicon} pb")

                    st.divider()
                    col1_res, col2_res = st.columns(2)
                    with col1_res:
                        st.subheader("Resultados Primer Forward")
                        st.metric("Posición (índice)", f"{posicion_fwd}")
                        st.metric("Distancia Total (Límite <= 5)", f"{dist_fwd}",
                                  delta=f"Pasa: {pasa_dist_total_fwd}", delta_color="off")
                        st.metric("Distancia Anclaje 3' (Límite 0)", f"{dist_3_fwd}",
                                  delta=f"Pasa: {pasa_anclaje_3_fwd}", delta_color="off")

                        st.markdown("**Alineamiento Forward:**")
                        label_primer = "Primer: "
                        label_target = "Target: "
                        padding = " " * len(label_primer)
                        alignment_str_fwd = generar_string_alineamiento(fwd_upper, subcadena_fwd)
                        display_string_fwd = f"{label_primer}{fwd_upper}\n{padding}{alignment_str_fwd}\n{label_target}{subcadena_fwd}"
                        st.code(display_string_fwd, language="text")

                    with col2_res:
                        st.subheader("Resultados Primer Reverse")
                        st.metric("Posición (índice)", f"{posicion_rev}")
                        st.metric("Distancia Total (Límite <= 5)", f"{dist_rev}",
                                  delta=f"Pasa: {pasa_dist_total_rev}", delta_color="off")
                        st.metric("Distancia Anclaje 3' (Límite 0)", f"{dist_3_rev}",
                                  delta=f"Pasa: {pasa_anclaje_3_rev}", delta_color="off")

                        st.markdown("**Primer (Input):**")
                        st.code(rev_upper, language="text")

                        st.markdown("**Alineamiento Reverse (RC vs Target):**")
                        label_rc = "RC:     "
                        label_target = "Target: "
                        padding = " " * len(label_rc)
                        alignment_str_rev = generar_string_alineamiento(rc_rev_upper, subcadena_rev)
                        display_string_rev = f"{label_rc}{rc_rev_upper}\n{padding}{alignment_str_rev}\n{label_target}{subcadena_rev}"
                        st.code(display_string_rev, language="text")

# =======================================================================
# --- PESTAÑA 2: GENERAR PRIMERS NUEVOS ---
# =======================================================================
with tab2:
    st.header("Generar Primers Nuevos (con Primer3)")

    st.subheader("Parámetros de Diseño del Producto")
    col1, col2 = st.columns(2)
    with col1:
        min_size = st.number_input(
            "Tamaño MÍNIMO del Producto (pb)",
            min_value=50, value=150, step=10, key="min_size_gen"
        )
    with col2:
        max_size = st.number_input(
            "Tamaño MÁXIMO del Producto (pb)",
            min_value=100, value=300, step=10, key="max_size_gen"
        )

    if st.button("¡Diseñar Ahora!", use_container_width=True, type="primary", key="btn_disenar"):
        rango_producto_usuario = [[int(min_size), int(max_size)]]
        st.header("Resultados del Diseño")
        st.info(f"Buscando 10 pares de primers con amplicones entre {min_size} y {max_size} pb...")

        try:
            seq_original = seq_input.strip().upper()
            seq_para_diseno = limpiar_degenerados(seq_original)

            if 'N' in seq_para_diseno and not 'N' in seq_original:
                st.warning(
                    "Advertencia: Se detectaron nucleótidos degenerados (V, S, R, etc.) en tu secuencia. Fueron convertidos a 'N' para el análisis con Primer3.")

            if len(seq_para_diseno) < 70:
                st.error(
                    f"Error: La secuencia ({len(seq_para_diseno)} pb) es más corta que el mínimo de 70 pb para el diseño.")
            else:
                resultados_primer3 = primer3.design_primers(
                    seq_args={'SEQUENCE_TEMPLATE': seq_para_diseno},
                    global_args={
                        'PRIMER_NUM_RETURN': 10,
                        'PRIMER_PRODUCT_SIZE_RANGE': rango_producto_usuario,
                        'PRIMER_OPT_SIZE': 20,
                        'PRIMER_MIN_SIZE': 18,
                        'PRIMER_MAX_SIZE': 22,
                        'PRIMER_OPT_TM': 60.0,
                        'PRIMER_MIN_TM': 58.0,
                        'PRIMER_MAX_TM': 62.0,
                        'PRIMER_MIN_GC': 40.0,
                        'PRIMER_MAX_GC': 60.0,
                        'PRIMER_EXPLAIN_FLAG': 1
                    }
                )
                num_pares_generados = resultados_primer3.get('PRIMER_PAIR_NUM_RETURNED', 0)
                if num_pares_generados > 0:
                    st.success(f"¡ÉXITO! Se generaron {num_pares_generados} pares de primers:")
                    amplicones_simulacion = np.zeros(num_pares_generados)
                    for i in range(num_pares_generados):

                        st.subheader(f"PAR DE PRIMERS # {i + 1}")
                        fwd_seq = resultados_primer3.get(f'PRIMER_LEFT_{i}_SEQUENCE')
                        rev_seq = resultados_primer3.get(f'PRIMER_RIGHT_{i}_SEQUENCE')
                        tamano_producto = resultados_primer3.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE')
                        amplicones_simulacion[i] = tamano_producto

                        try:
                            (GC_fwd_porc, GC_rev_porc, Tm_fwd, Tm_rev,
                             autocomp_fwd, autocomp_rev) = parametros_primer(fwd_seq, rev_seq)
                            col_fwd, col_rev = st.columns(2)
                            with col_fwd:
                                st.markdown("**Primer Forward**")
                                st.code(fwd_seq, language="text")
                                st.metric("Tm (Calculada)", f"{Tm_fwd:.2f} °C")
                                st.metric("%GC (Calculado)", f"{GC_fwd_porc:.2f} %")
                                st.metric("Autocomplementariedad (Dist. Levenshtein)", f"{autocomp_fwd}")
                                st.caption("Distancia: Un valor ALTO es bueno (poco palíndromo)")
                            with col_rev:
                                st.markdown("**Primer Reverse**")
                                st.code(rev_seq, language="text")
                                st.metric("Tm (Calculada)", f"{Tm_rev:.2f} °C")
                                st.metric("%GC (Calculado)", f"{GC_rev_porc:.2f} %")
                                st.metric("Autocomplementariedad (Dist. Levenshtein)", f"{autocomp_rev}")
                                st.caption("Distancia: Un valor ALTO es bueno (poco palíndromo)")
                            st.metric("Tamaño del Amplicón (de Primer3)", f"{tamano_producto} pb")
                            st.divider()
                        except Exception as e:
                            st.error(f"Error al calcular métricas personalizadas para el Par #{i + 1}: {e}")
                            st.divider()
                            continue

                    # Llamar a la simulación de electroforesis aquí para Pestaña 5
                    # Guardamos los amplicones para mostrarlos en la Pestaña 5
                    st.session_state["amplicones_simulacion_tab2"] = amplicones_simulacion.tolist()

                else:
                    st.error("FALLO EL DISEÑO")
                    st.text("Explicaciones de Primer3:")
                    if 'PRIMER_LEFT_EXPLAIN' in resultados_primer3:
                        st.text(f"  Errores (Forward): {resultados_primer3['PRIMER_LEFT_EXPLAIN']}")
                    if 'PRIMER_RIGHT_EXPLAIN' in resultados_primer3:
                        st.text(f"  Errores (Reverse): {resultados_primer3['PRIMER_RIGHT_EXPLAIN']}")
                    if 'PRIMER_PAIR_EXPLAIN' in resultados_primer3:
                        st.text(f"  Errores (Pares): {resultados_primer3['PRIMER_PAIR_EXPLAIN']}")
        except Exception as e:
            st.error(f"Ocurrió un error inesperado al ejecutar primer3: {e}")

# =======================================================================
# --- PESTAÑA 3: BÚSQUEDA DE MÚLTIPLES PARES DE PRIMERS ---
# =======================================================================
with tab3:
    st.header("Analizar Múltiples Pares de Primers (Batch)")
    st.info("Esta pestaña aplica el mismo filtro biológico: Distancia Total <= 5 Y Distancia Anclaje 3' == 0.")

    col1, col2 = st.columns(2)
    with col1:
        fwd_list_input = st.text_area(
            "Lista de Primers FORWARD (uno por línea)",
            value="ACAACAGCGTTGAAACAGCC\nAGCTTCTACAGGCTCCCAGT\n",
            height=250,
            key="fwd_list_input"
        )
    with col2:
        rev_list_input = st.text_area(
            "Lista de Primers REVERSE (uno por línea)",
            value="GGGGCTTTGAGCTTCGAGAT\nGGCTGTTTCAACGCTGTTGT\n",
            height=250,
            key="rev_list_input"
        )

    if st.button("Analizar Múltiples Pares", use_container_width=True, type="primary", key="btn_analizar_multi"):

        target_sequence_upper = seq_input.strip().upper()
        fwd_list_raw = fwd_list_input.strip().upper()
        rev_list_raw = rev_list_input.strip().upper()

        if not target_sequence_upper or not fwd_list_raw or not rev_list_raw:
            st.warning("Por favor, ingresa la secuencia target y ambas listas de primers.")
        elif len(target_sequence_upper) < 70:
            st.error(f"Error: El genoma de ({len(target_secuencia_upper)} pb) es más corto que el mínimo de 70 pb")
        else:
            fwd_list = [p.strip() for p in fwd_list_raw.split('\n') if p.strip()]
            rev_list = [p.strip() for p in rev_list_raw.split('\n') if p.strip()]

            if not fwd_list or not rev_list:
                st.error("Una de las listas de primers está vacía")
            elif len(fwd_list) != len(rev_list):
                st.error(
                    f"Error: El número de primers Forward ({len(fwd_list)}) no coincide con el número de primers Reverse ({len(rev_list)}).")
            else:
                st.header(f"Resultados del Análisis de {len(fwd_list)} Pares")
                amplicones_simulacion_batch = []
                with st.spinner("Analizando todos los pares..."):

                    num_exitos = 0
                    num_fallos = 0

                    for i in range(len(fwd_list)):
                        fwd_upper = fwd_list[i]
                        rev_upper = rev_list[i]

                        st.subheader(f"Par #{i + 1}: {fwd_upper[:10]}... / {rev_upper[:10]}...")

                        rc_rev_upper = reverso_complementario(rev_upper)
                        posicion_fwd, dist_fwd, subcadena_fwd = buscar_mejor_match_levenshtein(fwd_upper, target_sequence_upper)
                        posicion_rev, dist_rev, subcadena_rev = buscar_mejor_match_levenshtein(rc_rev_upper, target_sequence_upper)

                        if posicion_fwd == -1:
                            st.error("Error: El primer Forward no se encontró.")
                            num_fallos += 1
                        elif posicion_rev == -1:
                            st.error("Error: El Reverso Complementario del primer Reverse no se encontró.")
                            num_fallos += 1
                        elif posicion_fwd >= posicion_rev:
                            st.error(
                                f"Error: Fwd (pos {posicion_fwd}) se encontró DESPUÉS del Rev (pos {posicion_rev}).")
                            num_fallos += 1
                        else:
                            dist_3_fwd = distancia_levenshtein(fwd_upper[-5:], subcadena_fwd[-5:])
                            dist_3_rev = distancia_levenshtein(rc_rev_upper[-5:], subcadena_rev[-5:])

                            pasa_dist_total_fwd = (dist_fwd <= 5)
                            pasa_dist_total_rev = (dist_rev <= 5)
                            pasa_anclaje_3_fwd = (dist_3_fwd == 0)
                            pasa_anclaje_3_rev = (dist_3_rev == 0)

                            veredicto_final = (pasa_dist_total_fwd and pasa_dist_total_rev and
                                               pasa_anclaje_3_fwd and pasa_anclaje_3_rev)

                            if veredicto_final:
                                st.success("VEREDICTO: ¡ÉXITO! El par CUMPLE las reglas.")
                                num_exitos += 1
                                inicio_fwd = posicion_fwd
                                fin_rev = posicion_rev + len(rc_rev_upper)
                                tamano_amplicon = fin_rev - inicio_fwd
                                amplicones_simulacion_batch.append(tamano_amplicon)
                            else:
                                st.error("VEREDICTO: FALLO. El par NO cumple las reglas.")
                                num_fallos += 1
                                tamano_amplicon = "N/A (Fallo en Anclaje)"

                            st.metric("Tamaño Estimado del Amplicón", f"{tamano_amplicon} pb")

                            col1_res, col2_res = st.columns(2)
                            with col1_res:
                                st.metric("Posición Fwd", f"{posicion_fwd}")
                                st.metric("Distancia Total Fwd (<= 5)", f"{dist_fwd}")
                                st.metric("Distancia Anclaje 3' Fwd (== 0)", f"{dist_3_fwd}")

                                label_primer = "Primer: "
                                label_target = "Target: "
                                padding = " " * len(label_primer)
                                alignment_str_fwd = generar_string_alineamiento(fwd_upper, subcadena_fwd)
                                display_string_fwd = f"{label_primer}{fwd_upper}\n{padding}{alignment_str_fwd}\n{label_target}{subcadena_fwd}"
                                st.code(display_string_fwd, language="text")

                            with col2_res:
                                st.metric("Posición Rev", f"{posicion_rev}")
                                st.metric("Distancia Total Rev (<= 5)", f"{dist_rev}")
                                st.metric("Distancia Anclaje 3' Rev (== 0)", f"{dist_3_rev}")

                                label_rc = "RC:     "
                                label_target = "Target: "
                                padding = " " * len(label_rc)
                                alignment_str_rev = generar_string_alineamiento(rc_rev_upper, subcadena_rev)
                                display_string_rev = f"{label_rc}{rc_rev_upper}\n{padding}{alignment_str_rev}\n{label_target}{subcadena_rev}"
                                st.code(display_string_rev, language="text")

                        st.divider()

                # Guardamos los amplicones exitosos para la Pestaña 5
                if amplicones_simulacion_batch:
                    st.session_state["amplicones_simulacion_tab3"] = amplicones_simulacion_batch
                else:
                    st.session_state["amplicones_simulacion_tab3"] = []

# =======================================================================
# --- PESTAÑA 4: VISUALIZADOR MAPA GENÉTICO  ---
# =======================================================================
with tab4:
    st.header("🗺️ Visualizador Interactivo de Mapa Genético")
    st.markdown("Define los genes (CDS) y los segmentos amplificados (Amplicones/Primers) para visualizar el mapa.")

    # Longitud Total del Genoma
    longitud_total_mapa = st.number_input(
        "**Longitud Total del Genoma (bp):**",
        min_value=1,
        # Asegúrate de que este valor sea grande, o usa la longitud de tu secuencia:
        value=2969,  # Esto es lo más seguro si 'seq_input' tiene el genoma completo
        step=10,
        key="longitud_total_mapa"
    )
    if longitud_total_mapa <= 0:
        st.error("La longitud del genoma debe ser mayor a 0.")
        st.stop()

    st.subheader("📝 Entrada de CDS (Genes)")
    st.markdown("Define el **Nombre**, las coordenadas de **Inicio** y **Fin**, y la **Hebr**a (+1 o -1).")

    cds_default_data = [
        {"name": "NSs", "start": 69, "end": 1499, "strand": "+1"},
        {"name": "N", "start": 2193, "end": 2969, "strand": "-1"}
    ]

    cds_column_config = {
        "name": st.column_config.TextColumn("Nombre del CDS", required=True),
        "start": st.column_config.NumberColumn("Inicio (bp)", min_value=1, max_value=longitud_total_mapa,
                                               required=True),
        "end": st.column_config.NumberColumn("Fin (bp)", min_value=1, max_value=longitud_total_mapa, required=True),
        "strand": st.column_config.SelectboxColumn("Hebras", options=["+1", "-1"], required=True)
    }

    cds_data = st.data_editor(
        cds_default_data,
        column_config=cds_column_config,
        num_rows="dynamic",
        key="cds_editor_tab4"
    )

    st.subheader("🔬 Entrada de Amplicones (Primers)")
    st.markdown("Define el **Nombre** del amplicon y sus coordenadas de **Inicio** y **Fin**.")

    amplicon_default_data = [
        {"name": "Par 1 (76pb)", "start": 928, "end": 1003},
        {"name": "Par 2 (Test)", "start": 2300, "end": 2450}
    ]

    amplicon_column_config = {
        "name": st.column_config.TextColumn("Nombre del Amplicon", required=True),
        "start": st.column_config.NumberColumn("Inicio (bp)", min_value=1, max_value=longitud_total_mapa,
                                               required=True),
        "end": st.column_config.NumberColumn("Fin (bp)", min_value=1, max_value=longitud_total_mapa, required=True)
    }

    amplicon_data = st.data_editor(
        amplicon_default_data,
        column_config=amplicon_column_config,
        num_rows="dynamic",
        key="amplicon_editor_tab4"
    )

    st.subheader("Resultado del Mapa Genético")

    # Validaciones para el mapa
    datos_validos_mapa = True
    for entry in cds_data + amplicon_data:
        if entry.get("start") >= entry.get("end"):
            st.error(
                f"⚠️ Error: El valor de **Inicio** ({entry['start']}) debe ser menor al de **Fin** ({entry['end']}) en **{entry.get('name', 'una fila')}**.")
            datos_validos_mapa = False
            break

    if datos_validos_mapa and (cds_data or amplicon_data):
        try:
            # Generar el mapa
            figura = crear_mapa_genetico(longitud_total_mapa, cds_data, amplicon_data)
            st.pyplot(figura)

            # Opción para descargar la imagen
            buf = io.BytesIO()
            figura.savefig(buf, format="png", bbox_inches='tight')

            st.download_button(
                label="Descargar Mapa como PNG",
                data=buf.getvalue(),
                file_name="mapa_genetico.png",
                mime="image/png"
            )
        except Exception as e:
            st.error(f"Error al generar el mapa genético: {e}")

    elif not (cds_data or amplicon_data):
        st.warning("Agrega al menos un CDS o un Amplicon para generar el mapa.")

    else:
        st.info("Por favor, corrige los errores en los datos para visualizar el mapa.")

# =======================================================================
# --- PESTAÑA 5: SIMULAR ELECTROFORESIS  ---
# =======================================================================
with tab5:
    st.header("💻 Simulación de Electroforesis en Gel")
    st.markdown(
        "Visualiza la simulación de la electroforesis a partir de los **Amplicones generados o analizados** en las pestañas anteriores.")

    # Opción para cargar los amplicones de las pestañas
    opcion_amplicones = st.radio(
        "Fuente de Amplicones:",
        ("Carga Manual", "Usar Generados (Pestaña 2)", "Usar Analizados Exitosos (Pestaña 3)"),
        key="opcion_amplicones_tab5"
    )

    # Inicializar la lista de amplicones a usar
    amplicones_a_simular = []

    if opcion_amplicones == "Carga Manual":
        default_amplicones = "250, 500, 1500"
        amplicones_input_manual = st.text_area(
            "Tamaños de Amplicones (pb), separados por comas o nueva línea:",
            value=default_amplicones,
            key="amplicones_input_manual"
        )
        try:
            amplicones_str_list = [a.strip() for a in amplicones_input_manual.replace('\n', ',').split(',') if
                                   a.strip()]
            amplicones_a_simular = [int(a) for a in amplicones_str_list if a.isdigit()]
            if not amplicones_a_simular and amplicones_str_list:
                st.warning("Asegúrate de ingresar solo números enteros para los tamaños de amplicones.")
        except:
            st.error("Formato de entrada manual no válido.")
            amplicones_a_simular = []


    elif opcion_amplicones == "Usar Generados (Pestaña 2)":
        if "amplicones_simulacion_tab2" in st.session_state:
            amplicones_a_simular = st.session_state["amplicones_simulacion_tab2"]
            st.info(f"Cargados {len(amplicones_a_simular)} amplicones generados en la Pestaña 2.")
        else:
            st.warning("No hay amplicones generados en la Pestaña 2.")

    elif opcion_amplicones == "Usar Analizados Exitosos (Pestaña 3)":
        if "amplicones_simulacion_tab3" in st.session_state:
            amplicones_a_simular = st.session_state["amplicones_simulacion_tab3"]
            st.info(f"Cargados {len(amplicones_a_simular)} amplicones exitosos de la Pestaña 3.")
        else:
            st.warning("No hay amplicones exitosos analizados en la Pestaña 3.")

    # Parámetros adicionales para la simulación
    st.subheader("Parámetros de Simulación")
    col_conc, col_carriles = st.columns(2)
    with col_conc:
        concentracion_gel = st.slider(
            "Concentración de Agarosa (%)",
            min_value=0.5, max_value=2.0, value=1.0, step=0.1, key="concentracion_gel"
        )
    with col_carriles:
        num_carriles_max = 50  # Limitamos por visualización
        num_carriles = len(amplicones_a_simular) if amplicones_a_simular else 5
        st.metric("Número de Carriles (Automático)", num_carriles)

    if amplicones_a_simular:
        amplicones_np = np.array(amplicones_a_simular)

        st.subheader("Resultado de la Simulación")
        # Llamar a la función de simulación
        try:
            # La función simular_electroforesis necesita manejar el guardado y/o devolución de la imagen
            # Asumiendo que la función la guarda en un archivo o la devuelve. Aquí la ejecutamos:
            simular_electroforesis(amplicones_np, concentracion_gel)

            # Mostrar la imagen después de que la función la genere (asumiendo que guarda como jpg)
            # Nota: Necesitarás ajustar cómo simular_electroforesis maneja la salida de la imagen

            # NOTA: Tu código original tenía un bloque para mostrar la imagen guardada
            # Lo ajustamos ligeramente, pero depende de CÓMO esté implementada
            # la función 'simular_electroforesis' en tu archivo aparte.

            st.success("Simulación de electroforesis completada.")

            # --- Bloque de Visualización de Imagen (Si la función guarda la imagen) ---
            st.info(
                "La visualización aquí depende de cómo esté implementada la función `simular_electroforesis` en tu archivo aparte.")
            # Si 'simular_electroforesis' guarda la imagen en un archivo y necesitas mostrarla:
            # st.image(nombre_del_archivo_de_la_imagen, use_container_width=True)


        except Exception as e:
            st.error(f"Error al ejecutar la simulación de electroforesis: {e}")
            st.warning(
                "Asegúrate de que el archivo `electroforesis_gel.py` esté disponible y la función `simular_electroforesis` esté definida correctamente.")
    else:
        st.warning("No hay amplicones válidos para simular la electroforesis.")