import numpy as np
import streamlit as st
import primer3
from electroforesis_gel import simular_electroforesis
import matplotlib.pyplot as plt
import io
import primer_tools as prmt

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
tab1, tab2, tab3, tab4 = st.tabs([
    "🔬 Generar Primers Nuevos", # PESTAÑA 1
    "📊 Analizar Múltiples Pares de Primers", # PESTAÑA 2
    "🗺️ Visualizador Mapa Genético",  # PESTAÑA 3
    "💻 Simular Electroforesis"  # PESTAÑA 4
])

# =======================================================================
# --- PESTAÑA 1: GENERAR PRIMERS NUEVOS ---
# =======================================================================
with tab1:
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
            seq_para_diseno = prmt.limpiar_degenerados(seq_original)

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
                             autocomp_fwd, autocomp_rev) = prmt.parametros_primer(fwd_seq, rev_seq)
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

                    simular_electroforesis(amplicones_simulacion)

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
# --- PESTAÑA 2: BÚSQUEDA DE MÚLTIPLES PARES DE PRIMERS ---
# =======================================================================
with tab2:
    st.header("Analizar Múltiples Pares de Primers")
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

    if st.button("Analizar", use_container_width=True, type="primary", key="btn_analizar_multi"):

        target_sequence_upper = seq_input.strip().upper()
        fwd_list_raw = fwd_list_input.strip().upper()
        rev_list_raw = rev_list_input.strip().upper()

        if not target_sequence_upper or not fwd_list_raw or not rev_list_raw:
            st.warning("Por favor, ingresa la secuencia target y ambas listas de primers.")
        elif len(target_sequence_upper) < 70:
            st.error(f"Error: El genoma de ({len(prmt.target_secuencia_upper)} pb) es más corto que el mínimo de 70 pb")
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

                        rc_rev_upper = prmt.reverso_complementario(rev_upper)
                        posicion_fwd, dist_fwd, subcadena_fwd = prmt.buscar_mejor_match_levenshtein(fwd_upper, target_sequence_upper)
                        posicion_rev, dist_rev, subcadena_rev = prmt.buscar_mejor_match_levenshtein(rc_rev_upper, target_sequence_upper)

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
                            dist_3_fwd = prmt.distancia_levenshtein(fwd_upper[-5:], subcadena_fwd[-5:])
                            dist_3_rev = prmt.distancia_levenshtein(rc_rev_upper[-5:], subcadena_rev[-5:])

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
                                alignment_str_fwd = prmt.generar_string_alineamiento(fwd_upper, subcadena_fwd)
                                display_string_fwd = f"{label_primer}{fwd_upper}\n{padding}{alignment_str_fwd}\n{label_target}{subcadena_fwd}"
                                st.code(display_string_fwd, language="text")

                            with col2_res:
                                st.metric("Posición Rev", f"{posicion_rev}")
                                st.metric("Distancia Total Rev (<= 5)", f"{dist_rev}")
                                st.metric("Distancia Anclaje 3' Rev (== 0)", f"{dist_3_rev}")

                                label_rc = "RC:     "
                                label_target = "Target: "
                                padding = " " * len(label_rc)
                                alignment_str_rev = prmt.generar_string_alineamiento(rc_rev_upper, subcadena_rev)
                                display_string_rev = f"{label_rc}{rc_rev_upper}\n{padding}{alignment_str_rev}\n{label_target}{subcadena_rev}"
                                st.code(display_string_rev, language="text")

                        st.divider()

# =======================================================================
# --- PESTAÑA 3: VISUALIZADOR MAPA GENÉTICO  ---
# =======================================================================
with tab3:
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
    st.markdown("Define el **Nombre**, las coordenadas de **Inicio** y **Fin**, y la **Hebra** (+1 o -1).")

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
        {"name": "Par número 1", "start": 928, "end": 1003},
        {"name": "Par número 2", "start": 2300, "end": 2450}
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
            prmt.crear_mapa_genetico(longitud_total_mapa, cds_data, amplicon_data)
            st.image(f"./imagenes/mapa_genetico.jpg")

            with open(f"./imagenes/mapa_genetico.jpg", "rb") as file:
                st.download_button(
                    label="📥 Descargar Mapa Genético (.jpg)",
                    data=file,
                    file_name="mapa_genetico.jpg",
                    mime="image/jpeg"
                )
        except Exception as e:
            st.error(f"Error al generar el mapa genético: {e}")

    elif not (cds_data or amplicon_data):
        st.warning("Agrega al menos un CDS o un Amplicon para generar el mapa.")

    else:
        st.info("Por favor, corrige los errores en los datos para visualizar el mapa.")

# =======================================================================
# --- PESTAÑA 4: SIMULAR ELECTROFORESIS  ---
# =======================================================================
with tab4:
    st.subheader("Electroforesis generada")

    # Construir la ruta de la imagen
    output_file = f"./imagenes/Electroforesis.jpg"

    # Mostrar
    try:
        st.image(output_file, use_container_width=True)
    except:
        st.warning("Aún no se ha generado ninguna imagen de electroforesis.")
