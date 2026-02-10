Extractor y Convertidor de Coordenadas Químicas

Este programa automatiza la extracción de coordenadas cartesianas desde archivos de salida con extensión `.sumviz`, realiza la conversión de unidades de Bohr a Angstroms y genera archivos compatibles con software de visualización química.

Características

Extracción Inteligente: Localiza la sección "Nuclear Charges and Cartesian Coordinates" dentro del archivo de entrada.
Conversión de Unidades: Transforma automáticamente las coordenadas usando el factor de conversión 
1 Bohr = 0.529177 A˚
Generación de Formatos:
    XYZ: Archivo estándar para visualización de estructuras.
    SMILES: Genera la nomenclatura simplificada de la molécula utilizando la librería Open Babel.

Requisitos

Para ejecutar este script, necesitarás tener instalado:

1.  Python 3.x
2.  Open Babel: Librería fundamental para la interconversión de formatos químicos.
    bash
    En sistemas basados en Debian/Ubuntu
    sudo apt-get install python3-openbabel
    En sistemas Fedora/Nobara 
    sudo dnf install python3-openbabel

Estructura de Archivos

`archivo.sumviz`: Archivo de entrada con los datos de la simulación.
`resultado.smi`: Archivo de salida con la cadena SMILES generada.
`visualizacion.xyz`: Archivo de salida con las coordenadas en Angstroms.

Configuración

Antes de correr el programa, asegúrate de actualizar las rutas de las carpetas en el script:

python
carpeta = "/tu/ruta/de/trabajo/"
archivo_in = os.path.join(carpeta, "tu_archivo.sumviz")
