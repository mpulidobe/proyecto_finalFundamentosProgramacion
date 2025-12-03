import numpy as np
import streamlit as st
import primer3
from electroforesis_gel import simular_electroforesis
#holaaaaa
IUPAC_MAP = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'}, 'U': {'T'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'C', 'G'},
    'W': {'A', 'T'}, 'K': {'G', 'T'}, 'M': {'A', 'C'},
    'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'T'},
}

def compatibilidad(base1, base2):
    bases_set1 = IUPAC_MAP.get(base1, {base1})
    bases_set2 = IUPAC_MAP.get(base2, {base2})
    return bool(bases_set1.intersection(bases_set2))


def traduccion_ntdegenerados(sequence):
    validos = set('ACGTN')
    seq_limpia = []
    for base in sequence.upper():
        if base in validos:
            seq_limpia.append(base)
        else:
            seq_limpia.append('N')
    seq_para_diseno = "".join(seq_limpia)
    return seq_para_diseno


def reverso_complementario(sequence):
    """
    Calcula el reverso complementario de una secuencia de ADN,
    manejando correctamente los códigos degenerados IUPAC.
    """
    dicc_complementos = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
        'D': 'H', 'H': 'D', 'N': 'N',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
        'r': 'y', 'y': 'r', 's': 's', 'w': 'w',
        'k': 'm', 'm': 'k', 'b': 'v', 'v': 'b',
        'd': 'h', 'h': 'd', 'n': 'n'
    }
    sequence_invertida = sequence[::-1]
    list_complementos = [dicc_complementos.get(base, base) for base in sequence_invertida]
    return "".join(list_complementos)


def distancia_levenshtein(s1, s2):
    filas = len(s1) + 1
    columnas = len(s2) + 1
    matriz = [[0 for _ in range(columnas)] for _ in range(filas)]
    for i in range(filas):
        matriz[i][0] = i
    for j in range(columnas):
        matriz[0][j] = j
    for i in range(1, filas):
        for j in range(1, columnas):
            if compatibilidad(s1[i - 1], s2[j - 1]):
                costo_sustitucion = 0
            else:
                costo_sustitucion = 1
            eliminacion = matriz[i - 1][j] + 1
            insercion = matriz[i][j - 1] + 1
            sustitucion = matriz[i - 1][j - 1] + costo_sustitucion
            matriz[i][j] = min(eliminacion, insercion, sustitucion)
    return matriz[filas - 1][columnas - 1]


def smith_waterman_score(s1, s2, match=1, mismatch=-1, gap=-2):
    filas = len(s1) + 1
    cols = len(s2) + 1
    matrix = [[0 for _ in range(cols)] for _ in range(filas)]
    max_score = 0
    for i in range(1, filas):
        for j in range(1, cols):
            diagonal = matrix[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch)
            arriba = matrix[i - 1][j] + gap
            izquierda = matrix[i][j - 1] + gap
            current_score = max(0, diagonal, arriba, izquierda)
            matrix[i][j] = current_score
            if current_score > max_score:
                max_score = current_score
    return max_score


def match_levenshtein(primer, target_sequence):
    len_primer = len(primer)
    mejor_distancia = float('inf')
    mejor_match_subcadena = ""
    mejor_posicion = -1
    if len_primer > len(target_sequence):
        return -1, float('inf'), ""
    ancho_ventana = len_primer
    for i in range(len(target_sequence) - len_primer + 1):
        subcadena = target_sequence[i: i + ancho_ventana]
        distancia = distancia_levenshtein(primer, subcadena)
        if distancia < mejor_distancia:
            mejor_distancia = distancia
            mejor_match_subcadena = subcadena
            mejor_posicion = i
            if distancia == 0:
                break
    return mejor_posicion, mejor_distancia, mejor_match_subcadena


def generar_string_alineamiento(s1, s2):
    alineamiento = []
    for i in range(len(s1)):
        # Usamos la misma lógica "inteligente" de compatibilidad
        if compatibilidad(s1[i], s2[i]):
            alineamiento.append('|')
        else:
            alineamiento.append(' ')
    return "".join(alineamiento)


def GC_Tm_autocomplemetary(primer_fwd, primer_rev):
    G_fwd_count = primer_fwd.count('G')
    C_fwd_count = primer_fwd.count('C')
    GC_fwd = G_fwd_count + C_fwd_count
    GC_fwd_porc = (GC_fwd / len(primer_fwd)) * 100
    G_rev_count = primer_rev.count('G')
    C_rev_count = primer_rev.count('C')
    GC_rev = G_rev_count + C_rev_count
    GC_rev_porc = (GC_rev / len(primer_rev)) * 100
    Tm_fwd = 64.9 + (41 * (GC_fwd - 16.4)) / len(primer_fwd)
    Tm_rev = 64.9 + (41 * (GC_rev - 16.4)) / len(primer_rev)
    rc_fwd = reverso_complementario(primer_fwd)
    autocomplemetary_fwd = distancia_levenshtein(primer_fwd, rc_fwd)
    rc_rev = reverso_complementario(primer_rev)
    autocomplemetary_rev = distancia_levenshtein(primer_rev, rc_rev)
    return round(GC_fwd_porc, 2), round(GC_rev_porc, 2), round(Tm_fwd, 2), round(Tm_rev,
                                                                                 2), autocomplemetary_fwd, autocomplemetary_rev


# --- INTERFAZ DE STREAMLIT ---

st.set_page_config(page_title="Analizador de Primers", layout="wide")
st.title("Primer4")

# --- ENTRADA DE SECUENCIA PRINCIPAL ---
st.header("Entrada de Secuencia Target")
default_sequence = """GTCATCTTTTATTCTTAATCAAACCTCACTCAGAAAAATCCAGAACTGTNNNAATAAGAACCAAAGCCACAVSRATGTCTACAAGTGTTTATGAATCGATCATTCAGACGAAAGCTTCGGTCTGGGGATCTACTGCATCTGGCAAATCTATTGTGGACTCTTACTGGATCCACGAGTTTTTAACTGGTTCTCCATTGATTCAAACTCAGTTGTATTCTGATTCAAGAAGCAAAAGCAGCTTTGGCTACACCACACGAGTTGGTGATCTTCVDSBCTTCAGAAGAGAAAGAGATTCTTTCTCAACACTTGTACATCCCTATTTTTGATGACATTGATTTCAACATCAATATCAATGATTCAGTCATGACAGTATCCGTTTGCTCCAACACGGTCAATGCTAATGGAGTGAAACATCAGGGTCATCTGAAGGTGCTTTCTCTTGCTCAACTGCACTCTATAGAGCCTACAATGAGCCGATCTGACATTGCTGACAGATTCCGTCTTCAAGAAAAAGACGTGATTCCCAATGACAGATACATTGATGCTGCTAACAAAGGCTCTCTTTCATGTGTTAAAGAGCATTCCTATAAAGTCGAAATGTGCCACAACCAAGCATTAGGTAAAGTTAATGTTCTATCCCCTAACAGAAATGTTCATGAATGGCTGTACAGCTTCAAGCCAGCTTTCAACCAGATTGAAAGCAACAACAGAACTGTAAATTCTCTTGCAGTGAAATCTCTGCTCATGTCTGCAGAAAACAACATAATGCCTAACTCTCAGGCCTTTGTCAAAGCTTCTACAGGCTCCCAGTTCAAGCTAAACCTCTGGCTGAGGATTCCTAAAGTTCTGAAACAGGTTTCTATTCAGAAACTATTTAAAGTTGCAGAAGATGAAACAGACAAAAGCTTTTATTTGTCTATTGCTTGCATCCCTAACCACAACAGCGTTGAAACAGCVSTMCTTGAATGTGACCATCATCTGCAAGCATCAGCTCCCAATCTCGAAGCTCAAAGCCCCTTTWKDTGAATTGACAATGATGTTTTCTGATCTGAGAGAACCTTACAACGTTGTGCATGATCCTTCTTACCCTCAAAGAATTGTTCATGCTCTGCTTGAGACACACACATCTTTTGCCCAAACTCTTTGCAATAACTTGCAAGAAGATGTGGTCATCTACACTTTGAACAACCCTGAGCTGACTTCTTTAAAGTTAGATTTAGGTAAGAAAACCCTAAATTACAGTGAAGATGCTTATAATAAGAAATATTTTCTTTCAAAAACTCTTGAATGCCTCCCAGTAAACACACAGACTATGTCTTATTTAGACAGCATTCAAATTCCCTCATGGAAGATTGACTTTGCCAGAGGAGAAATCAAAATTTCCCCTCAATCAATCTCTGTTGCAAAATCTTTgttgaagctagatcttgatgtgatcagaggaaagaaatctctgcctcagggagcttctgaatcagagtcaaagcaatttgtgtctatttgtctgctcctttaactattcctttcctctttaaatcttctttcagcttctttccaacttctttcaatctcttttcargtttctctcatttcctttaaatttcctaattttctttacttcctttct"""

seq_input = st.text_area("Secuencia de Referencia (Target)",
                         value=default_sequence,
                         height=250, key="seq_input_main")

st.divider()

# --- NAVEGACIÓN POR PESTAÑAS ---
tab1, tab2, tab3, tab4 = st.tabs([
    "🧬 Analizar Par de Primers (PCR)",
    "🔬 Generar Primers Nuevos",
    "📊 Analizar Múltiples Pares (Batch)",
    "💻 Simular Electroforesis"
])

# =======================================================================
# --- PESTAÑA 1: ANALIZAR PAR ---
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
                    posicion_fwd, dist_fwd, subcadena_fwd = match_levenshtein(fwd_upper, seq_upper)
                    posicion_rev, dist_rev, subcadena_rev = match_levenshtein(rc_rev_upper, seq_upper)

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

                        ### ¡VISUALIZACIÓN CORREGIDA! ###
                        st.markdown("**Alineamiento Forward:**")
                        # Definimos etiquetas con padding
                        label_primer = "Primer: "  # 8 chars
                        label_target = "Target: "  # 8 chars
                        padding = " " * len(label_primer)  # 8 espacios
                        # Generamos el string de |
                        alignment_str_fwd = generar_string_alineamiento(fwd_upper, subcadena_fwd)
                        # Creamos el string final limpio
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

                        ### ¡VISUALIZACIÓN CORREGIDA! ###
                        st.markdown("**Alineamiento Reverse (RC vs Target):**")
                        # Definimos etiquetas con padding
                        label_rc = "RC:     "  # 8 chars
                        label_target = "Target: "  # 8 chars
                        padding = " " * len(label_rc)  # 8 espacios
                        # Generamos el string de |
                        alignment_str_rev = generar_string_alineamiento(rc_rev_upper, subcadena_rev)
                        # Creamos el string final limpio
                        display_string_rev = f"{label_rc}{rc_rev_upper}\n{padding}{alignment_str_rev}\n{label_target}{subcadena_rev}"
                        st.code(display_string_rev, language="text")

# =======================================================================
# --- PESTAÑA 2: GENERAR PRIMERS ---
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
            seq_para_diseno = traduccion_ntdegenerados(seq_original)

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
                             autocomp_fwd, autocomp_rev) = GC_Tm_autocomplemetary(fwd_seq, rev_seq)
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
                    simular_electroforesis(amplicones_simulacion, 1)
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
# --- PESTAÑA 3: BÚSQUEDA DE MÚLTIPLES PARES ---
# =======================================================================
with tab3:
    st.header("Analizar Múltiples Pares de Primers (Batch)")
    st.info("Esta pestaña aplica el mismo filtro biológico: Distancia Total <= 5 Y Distancia Anclaje 3' == 0.")

    col1, col2 = st.columns(2)
    with col1:
        fwd_list_input = st.text_area(
            "Lista de Primers FORWARD (uno por línea)",
            value="ACAACAGCGTTGAAACAGCC\nAGCTTCTACAGGCTCCCAGT\nPRIMER_FWD_3",
            height=250,
            key="fwd_list_input"
        )
    with col2:
        rev_list_input = st.text_area(
            "Lista de Primers REVERSE (uno por línea)",
            value="GGGGCTTTGAGCTTCGAGAT\nGGCTGTTTCAACGCTGTTGT\nPRIMER_REV_3",
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
            st.error(f"Error: El genoma de ({len(target_sequence_upper)} pb) es más corto que el mínimo de 70 pb")
        else:
            fwd_list = [p.strip() for p in fwd_list_raw.split('\n') if p.strip()]
            rev_list = [p.strip() for p in rev_list_raw.split('\n') if p.strip()]

            if not fwd_list or not rev_list:
                st.error("Una de las listas de primers está vacía.")
            elif len(fwd_list) != len(rev_list):
                st.error(
                    f"Error: El número de primers Forward ({len(fwd_list)}) no coincide con el número de primers Reverse ({len(rev_list)}).")
            else:
                st.header(f"Resultados del Análisis de {len(fwd_list)} Pares")
                with st.spinner("Analizando todos los pares..."):

                    num_exitos = 0
                    num_fallos = 0

                    for i in range(len(fwd_list)):
                        fwd_upper = fwd_list[i]
                        rev_upper = rev_list[i]

                        st.subheader(f"Par #{i + 1}: {fwd_upper[:10]}... / {rev_upper[:10]}...")

                        rc_rev_upper = reverso_complementario(rev_upper)
                        posicion_fwd, dist_fwd, subcadena_fwd = match_levenshtein(fwd_upper, target_sequence_upper)
                        posicion_rev, dist_rev, subcadena_rev = match_levenshtein(rc_rev_upper, target_sequence_upper)

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
                            else:
                                st.error("VEREDICTO: FALLO. El par NO cumple las reglas.")
                                num_fallos += 1

                            inicio_fwd = posicion_fwd
                            fin_rev = posicion_rev + len(rc_rev_upper)
                            tamano_amplicon = fin_rev - inicio_fwd
                            st.metric("Tamaño Estimado del Amplicón", f"{tamano_amplicon} pb")

                            col1_res, col2_res = st.columns(2)
                            with col1_res:
                                st.metric("Posición Fwd", f"{posicion_fwd}")
                                st.metric("Distancia Total Fwd (<= 5)", f"{dist_fwd}")
                                st.metric("Distancia Anclaje 3' Fwd (== 0)", f"{dist_3_fwd}")

                                ### ¡VISUALIZACIÓN CORREGIDA! ###
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

                                ### ¡VISUALIZACIÓN CORREGIDA! ###
                                label_rc = "RC:     "
                                label_target = "Target: "
                                padding = " " * len(label_rc)
                                alignment_str_rev = generar_string_alineamiento(rc_rev_upper, subcadena_rev)
                                display_string_rev = f"{label_rc}{rc_rev_upper}\n{padding}{alignment_str_rev}\n{label_target}{subcadena_rev}"
                                st.code(display_string_rev, language="text")

                        st.divider()

# =======================================================================
# --- PESTAÑA 4: VISUALIZACIÓN DE ELECTROFORESIS ---
# =======================================================================
with tab4:
    st.subheader("Electroforesis generada")

    # Construir la ruta de la imagen
    numero_imagen = st.session_state.get("numero_imagen", 1)
    output_file = f"../imagen/Electroforesis{numero_imagen}.jpg"

    # Mostrar
    try:
        st.image(output_file, use_container_width=True)
    except:
        st.warning("Aún no se ha generado ninguna imagen de electroforesis.")
