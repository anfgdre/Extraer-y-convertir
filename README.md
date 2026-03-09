# **Extraer y Convertir**
Autores:
- [Andrea](anfgdre)

Este proyecto automatiza la extracción de coordenadas cartesianas desde archivos `.sumviz` (en unidades Bohr) y las convierte a formato XYZ (con unidades Angstrom). Adicionalmente, esta librería obtiene el identificador químico `SMILES` utilizando la librería Open Babel. 

Este un proyecto desarrollado para el servicio social en el área de quimioinformática.

## Estructura del Proyecto

* `src/`: Contiene el script principal de procesamiento (`convertidor.py`).
* `data/input/`: Carpeta donde se deben colocar los archivos `.sumviz` originales.
* `data/output/`: Carpeta donde se guardarán los resultados (`.smi` y `.xyz`).
* `tests/`: Scripts de validación de datos de entrada.

## Requisitos e instalación
Para el uso de esta librería se sugiere la creación del ambiente conda correspondiente.

Python 3.x
Open Babel: Necesario para la generación de SMILES y archivos XYZ.

## Instrucciones de Uso

1. Preparar los datos
Coloca tus archivos con extensión `.sumviz` en la carpeta `data/input/`.

2. Validar los archivos
Antes de procesar, verifica que los archivos tengan el formato correcto ejecutando el script de prueba de la siguiente manera:

```bash
python tests/test_inputs.py
```