import os

def validar_archivos_sumviz():
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    INPUT_DIR = os.path.join(BASE_DIR, "data", "input")
    
    archivos = [f for f in os.listdir(INPUT_DIR) if f.endswith(".sumviz")]
    
    if not archivos:
        print("No se encontraron archivos .sumviz para procesar.")
        return

    for archivo in archivos:
        ruta = os.path.join(INPUT_DIR, archivo)
        conteo_atomos = 0
        seccion_encontrada = False
        
        with open(ruta, 'r') as f:
            contenido = f.read()
            if "Nuclear Charges and Cartesian Coordinates" in contenido:
                seccion_encontrada = True
        
        if seccion_encontrada:
            print(f"{archivo}: Formato válido.")
        else:
            print(f"{archivo}: Falta la sección de coordenadas.")

if __name__ == "__main__":
    validar_archivos_sumviz()