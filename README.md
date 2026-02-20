Extraer y Convertir 

Este proyecto automatiza la extracción de coordenadas cartesianas desde archivos `.sumviz` (formato Bohr) y las convierte a formato XYZ (Angstrom) y SMILES utilizando la librería Open Babel. Es un proyecto desarrollado para el servicio social en el área de quimioinformática.

Estructura del Proyecto

* `src/`: Contiene el script principal de procesamiento (`convertidor.py`).
* `data/input/`: Carpeta donde se deben colocar los archivos `.sumviz` originales.
* `data/output/`: Carpeta donde se guardarán los resultados (`.smi` y `.xyz`).
* `tests/`: Scripts de validación de datos de entrada.

equisitos

Python 3.x
Open Babel: Necesario para la generación de SMILES y archivos XYZ.

Instrucciones de Uso

1. Preparar los datos
Coloca tus archivos con extensión `.sumviz` en la carpeta `data/input/`.

2. Validar los archivos
Antes de procesar, verifica que los archivos tengan el formato correcto ejecutando el script de prueba:
```bash
python tests/test_inputs.py