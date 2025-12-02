import cv2
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')



def simular_electroforesis(amplicones_array, numero_imagen):
    ''' Tipos de matrices (arrays) de numpy:
    Array unidimensional 1D (N) 'son como los ejes' en este caso las distancias
    Array bidimensional 2D (N, M) 'el espacio/superficie de la imagen' en este caso el gel

    Función meshgrid
    Permite generar dos arrays bidimensionales a partir de dos arrays unidimensionales.
    A partir de las matrices unidimensionales x y y, genero los siguientes arrays bidimensionales que llamé 'malla_':
    malla_x[f, c] posición horizontal
    malla_y[c, f] posición vertical'''

    '''Conversión del array de coordenadas de pixeles a un array de distancias físicas
    Para ello utilizo las dimensiones reales del gel y del marcador comercial utilizados:
    El largo del gel es de 25 cm y el marcador contiene fragmentos entre 200-6200 pb con incrementos de 200pb
    El valor mínimo y máximo de las coordenadas de pixeles de los margenes de dicho marcador, 
    corresponden a 0 y 250 mm en distancia física respectivamente. 
    La distancia física de los pixeles ubicados entre esos valores en el eje y, se determinan utilizando una interpolación lineal 
    Se crearon dos array de numpy uno con los valores (en pb) y otro con su correspondiente distancia física en (mm)'''

    marcador_bp = np.array([1500, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100])

    distancia_marcador = np.array([27, 50, 57, 65, 74, 83, 95, 109, 126, 147, 175])

    '''Para determinar la distancia en mm de una muestra se utiliza un modelo logaritmico 
    (pag 255 del libro Principles of Physical Biochemistry)'''
    log_bp = np.log10(marcador_bp)
    coeficientes_lineaTendencia = np.polyfit(log_bp, distancia_marcador, 7)


    '''Dimensiones físicas del gel'''
    ancho_gel = 200      # [mm]
    alto_gel = 250       # [mm]


    '''Malla de coordenadas (en mm)
    Se tomaron 1500 puntos como la cantidad de puntos en los que se dividen la region tanto en el eje x como en y'''
    puntos = 1500
    x = np.linspace(0, ancho_gel, puntos)
    y = np.linspace(0, alto_gel, puntos)

    '''Se utiliza la función np.meshgrid para crear una malla bidimensional a partir de dos arrays unidimensionales en x y y'''
    malla_x, malla_y = np.meshgrid(x, y)


    '''Se prepara el array de intensidades 2D (cuyos valores pueden estar entre 0 y 1) y va a ser 1 cuando encuentre una banda 
    para ello se crea un array de ceros con la misma cantidad de posiciones del meshgrid que luego será reemplazado, se puede usar malla_x o malla_y
    como entrada porque ambos tienen el mismo shape'''
    intensidades_array = np.zeros_like(malla_x)


    '''Se utiliza un peine que contiene 16 pocillos y el anchos de cada pocillo corresponde a 1.5 mm
    El ancho del gel es de 20 cm pero se respeta una margen de 10 mm en los bordes laterales'''
    cantidad_pocillos = 16
    ancho_pocillo = 1.5
    marcador_pocillo = 6
    margen_gel = 10


    '''Teniendo en cuenta esas margenes, se calcula el ancho disponible y se encuentra la distancia entre un pocillo y otro
    dividiendo ese ancho disponible entre el número de pocillos'''
    ancho_disponible = ancho_gel - 2*margen_gel
    distancia_entrePocillos = ancho_disponible / (cantidad_pocillos - 1) #Distancia entre un pocillo y otro


    '''Se prepara el array 1D en el que se van a guardar los valores de las coordenadas en x de cada pocillo, 
    para ello de crea un array de ceros que luego va a ser reemplazado con dichos valores'''
    centro_pocillo = np.zeros(cantidad_pocillos)

    '''Teniendo en cuenta que son 16 pocillos 
    Se calcular la coordenada en y para cada pocillo, teniendo en cuenta que '''
    for pocillo in range(cantidad_pocillos):
        centro_pocillo[pocillo] = margen_gel + pocillo * distancia_entrePocillos


    '''Definir funciones para añadir una banda'''

    def agregar_banda(distancia_banda, ubicacionBanda_x, grosor_banda=1):
        """Agrega una banda rectangular basada en distancias reales en mm."""
        condicion_vertical   = np.abs(malla_y - distancia_banda) < (grosor_banda/2) #cuando esta condición se cumple da true
        condicion_horizontal = np.abs(malla_x - ubicacionBanda_x ) < (ancho_pocillo/2) #cuando esta condición se cumple da true
        intensidades_array[np.logical_and(condicion_vertical, condicion_horizontal)] = 1.0

    def dibujar_pocillo(ubicacionPocillo_x):
        """Agrega un pocillo físico."""
        condicion_y = (malla_y < marcador_pocillo)
        condicion_x = np.abs(malla_x - ubicacionPocillo_x) < (ancho_pocillo / 2)
        intensidades_array[np.logical_and(condicion_y, condicion_x)] = 0.5   # intensidad media (color gris)

    '''Dibujar pocillos (cuadrados grises)'''
    for i in centro_pocillo: #i va de 0 a 16
        dibujar_pocillo(i)

    '''El marcador está ubicado en el primer pocillo (dibuje bandas blancas en el primer carril)'''
    pocillo_marcador = 0
    x_marc = centro_pocillo[pocillo_marcador]

    for d in distancia_marcador:
        agregar_banda(d, x_marc)

    for counter in range(len(amplicones_array)):
        amplicon_bp = amplicones_array[counter]
        distancia_amplicon = np.polyval(coeficientes_lineaTendencia, np.log10(amplicon_bp))

        '''Las muestras ubicadas en los carriles de posiciones >= 1'''
        pocillo_muestra = counter + 1
        x_sample = centro_pocillo[pocillo_muestra]
        agregar_banda(distancia_amplicon, x_sample)

    imagen_blurr = cv2.GaussianBlur(intensidades_array, (5, 5), 0)

    plt.figure(figsize=(13, 9))
    plt.imshow(imagen_blurr, cmap="gray",
               extent=[0, ancho_gel, alto_gel, 0], aspect='auto')

    ax = plt.gca()  # eje principal

    # ===========================
    # 1. Mover eje X arriba
    # ===========================
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()

    # ===========================
    # 2. Eje Y izquierdo: tamaño en bp
    # ===========================
    ax.set_yticks(distancia_marcador)
    ax.set_yticklabels([f"{bp}" for bp in marcador_bp])
    plt.ylabel("Tamaño (bp)")

    # ===========================
    # 3. Eje X arriba: etiquetas de pocillos
    # ===========================
    etiquetas_x = ["m"] + [str(i) for i in range(1, cantidad_pocillos)]
    ax.set_xticks(centro_pocillo)
    ax.set_xticklabels(etiquetas_x, rotation=0)

    # ===========================
    # 4. Eje derecho: distancia mm
    # ===========================
    ax2 = ax.twinx()  # segundo eje Y
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks(distancia_marcador)
    ax2.set_yticklabels([f"{d:.1f}" for d in distancia_marcador])
    ax2.set_ylabel("Migración (mm)")

    plt.tight_layout()
    output_file = f'../imagen/Electroforesis{numero_imagen}.jpg'
    plt.savefig(output_file)

