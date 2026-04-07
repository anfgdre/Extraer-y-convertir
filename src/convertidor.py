import os
import glob
import math
import time
import json
from collections import Counter
from datetime import datetime
# Asegúrate de tener instalado openbabel (pip install openbabel)
from openbabel import openbabel

# --- CONFIGURACIÓN DE RUTAS ---
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUT_DIR = os.path.join(BASE_DIR, "data", "input")
OUTPUT_BASE_DIR = os.path.join(BASE_DIR, "data", "output")

def calcular_distancia(p1, p2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))

def obtener_metadatos_molecula(lista_coords):
    puntos = []
    atomos = []
    for c in lista_coords:
        parts = c.split()
        atomos.append(parts[0])
        puntos.append(tuple(map(float, parts[1:])))
    
    max_d = 0.0
    for i in range(len(puntos)):
        for j in range(i + 1, len(puntos)):
            d = calcular_distancia(puntos[i], puntos[j])
            if d > max_d: max_d = d
            
    formula = dict(Counter(atomos))
    return round(max_d, 3), formula

def extraer_coordenadas_sumviz(ruta_archivo):
    coordenadas = []
    en_seccion = False
    factor_bohr_to_angstrom = 0.529177
    elementos_quimicos = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", 
                         "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca"}

    try:
        with open(ruta_archivo, 'r', encoding='utf-8') as f:
            for linea in f:
                linea_limpia = linea.strip()
                if "Nuclear Charges and Cartesian Coordinates" in linea_limpia:
                    en_seccion = True
                    continue
                if en_seccion:
                    if not linea_limpia and len(coordenadas) > 0:
                        en_seccion = False
                        break
                    if "Some Atomic Properties" in linea_limpia or "q(A) =" in linea_limpia:
                        en_seccion = False
                        break
                    partes = linea_limpia.split()
                    if len(partes) == 5:
                        nombre_raw = partes[0]
                        simbolo = "".join([c for c in nombre_raw if c.isalpha()])
                        if simbolo in elementos_quimicos:
                            try:
                                x = float(partes[2]) * factor_bohr_to_angstrom
                                y = float(partes[3]) * factor_bohr_to_angstrom
                                z = float(partes[4]) * factor_bohr_to_angstrom
                                coordenadas.append(f"{simbolo} {x:10.6f} {y:10.6f} {z:10.6f}")
                            except ValueError:
                                continue
        if not coordenadas: return None, 0, {}
        max_dist, formula = obtener_metadatos_molecula(coordenadas)
        return coordenadas, max_dist, formula
    except Exception:
        return None, 0, {}

def ejecutar_procesamiento(incluir_smiles=False):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    nombre_carpeta = f"Ejecucion_{timestamp}"
    ruta_salida_actual = os.path.join(OUTPUT_BASE_DIR, nombre_carpeta)
    os.makedirs(ruta_salida_actual, exist_ok=True)

    archivos = glob.glob(os.path.join(INPUT_DIR, "*.sumviz"))
    if not archivos:
        print(f"\n[!] No se encontraron archivos en: {INPUT_DIR}")
        return

    reporte_maestro = {
        "sesion": nombre_carpeta,
        "fecha_ejecucion": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "total_archivos": len(archivos),
        "resultados": []
    }

    for ruta in archivos:
        start_time = time.time()
        nombre = os.path.splitext(os.path.basename(ruta))[0]
        coord_list, max_d, formula = extraer_coordenadas_sumviz(ruta)
        
        datos = {"archivo": nombre, "estado": "Error", "metadatos": {}, "conversion": {}}

        if coord_list:
            xyz_block = f"{len(coord_list)}\n{nombre}\n" + "\n".join(coord_list)
            with open(os.path.join(ruta_salida_actual, f"{nombre}.xyz"), "w") as f:
                f.write(xyz_block)
            
            datos["metadatos"] = {"distancia_max_A": max_d, "conteo_atomos": formula}
            
            if incluir_smiles:
                obConv = openbabel.OBConversion()
                obConv.SetInAndOutFormats("xyz", "smi")
                mol = openbabel.OBMol()
                if obConv.ReadString(mol, xyz_block):
                    mol.ConnectTheDots()
                    mol.PerceiveBondOrders()
                    smiles = obConv.WriteString(mol).strip()
                    with open(os.path.join(ruta_salida_actual, f"{nombre}.smi"), "w") as f:
                        f.write(smiles + "\n")
                    datos["estado"] = "Exito"
                    datos["conversion"]["smiles"] = smiles
                else:
                    datos["estado"] = "Error SMILES"
            else:
                datos["estado"] = "Exito"
        
        datos["conversion"]["tiempo_seg"] = round(time.time() - start_time, 4)
        reporte_maestro["resultados"].append(datos)
        print(f"-> {nombre}: {datos['estado']}")

    with open(os.path.join(ruta_salida_actual, "reporte.json"), 'w') as f:
        json.dump(reporte_maestro, f, indent=4)

def menu():
    while True:
        print("\n1. Extraer XYZ\n2. Extraer XYZ + SMILES\n3. Salir")
        op = input("Opción: ")
        if op == "1": ejecutar_procesamiento(False)
        elif op == "2": ejecutar_procesamiento(True)
        elif op == "3": break

if __name__ == "__main__":
    menu()