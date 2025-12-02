import cv2 as cv
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.use('TkAgg')

'''En análisis de imagen, los valores de intensidad (rango dinámico) van de 0 a 255
Son 256 niveles
El valor 0 representa el color negro
EL valor 255 representa el color blanco
La función lectura_imagen lee el archivo de imagen en escala de grises y lo convierte en un array de pixeles
Esto se hizo porque la entrada de la función plt.imshow (Matplotlib) es un array, no un archivo .jpg
Matplotlib identifica el valor de intensidad de cada pixel (eje z) y con esta información puedo obtener las coordenadas de cada pixel'''

def lectura_imagen(path):
    gray_array = cv.imread(path, cv.IMREAD_GRAYSCALE)
    return gray_array

def graficar_imagen(imagen):
    plt.imshow(imagen, cmap='gray')
    plt.show()

imagen_array = lectura_imagen("../imagen/marcador_pb.jpg")
graficar_imagen(imagen_array)