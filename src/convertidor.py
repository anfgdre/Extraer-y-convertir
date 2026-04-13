import os
import glob
import json
import math
import shutil
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from openbabel import openbabel
from datetime import datetime
from collections import Counter #

# --- CONFIGURACIÓN DE RUTAS ---
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUT_DIR = os.path.join(BASE_DIR, "data", "input")
OUTPUT_BASE_DIR = os.path.join(BASE_DIR, "data", "output")

def generar_readme(ruta_carpeta):
    """Crea la documentación técnica del proceso."""
    contenido = f"""
SISTEMA DE PROCESAMIENTO MOLECULAR + AUDITORÍA (RDKit/OpenBabel)
Fecha: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
--------------------------------------------------
- Filtro de Calidad: Moléculas con anomalías o fallos se mueven a /debug.
- Validación Cruzada: Se comparan SMILES de OpenBabel y RDKit.
- Test Visual: Se genera un histograma de distancias por cada archivo.
- Detección de Anomalías: Identifica átomos vecinos a menos de 0.8 Å.
- Conteo de Átomos: Clasificación por tipo y cantidad.
--------------------------------------------------
"""
    with open(os.path.join(ruta_carpeta, "README_PROCESO.txt"), "w", encoding="utf-8") as f:
        f.write(contenido)

def extraer_datos(ruta_archivo):
    """Extrae coordenadas y convierte de Bohr a Ángstrom."""
    coordenadas, puntos, atomos = [], [], []
    factor = 0.529177
    en_seccion = False
    elementos = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca"}

    try:
        with open(ruta_archivo, 'r', encoding='utf-8') as f:
            for linea in f:
                if "Nuclear Charges and Cartesian Coordinates" in linea:
                    en_seccion = True; continue
                if en_seccion:
                    partes = linea.split()
                    if len(partes) == 5:
                        simbolo = "".join([c for c in partes[0] if c.isalpha()])
                        if simbolo in elementos:
                            x, y, z = float(partes[2])*factor, float(partes[3])*factor, float(partes[4])*factor
                            coordenadas.append(f"{simbolo} {x:10.6f} {y:10.6f} {z:10.6f}")
                            puntos.append((x, y, z))
                            atomos.append(simbolo)
                    elif len(coordenadas) > 0: break
        return coordenadas, puntos, atomos
    except: return None, None, None

def calcular_matriz_y_max(puntos, atomos):
    """Calcula distancias e identifica anomalías por cercanía extrema (< 0.8 Å)."""
    matriz = []
    lista_dist = []
    anomalias_vecinos = 0 #
    n = len(puntos)
    for i in range(n):
        for j in range(i + 1, n):
            d = math.sqrt(sum((a - b) ** 2 for a, b in zip(puntos[i], puntos[j])))
            
            # Conteo de anomalías: átomos demasiado cerca (choque atómico)
            if d < 0.8:
                anomalias_vecinos += 1
            
            matriz.append({"de": f"{atomos[i]}({i})", "a": f"{atomos[j]}({j})", "distancia_A": round(d, 4)})
            lista_dist.append(d)
            
    dist_max = max(lista_dist) if lista_dist else 0
    return matriz, lista_dist, dist_max, anomalias_vecinos

def generar_histograma(distancias, nombre, ruta_carpeta):
    """Genera test visual de distancias."""
    if not distancias: return
    plt.figure(figsize=(8, 5))
    plt.hist(distancias, bins=15, color='skyblue', edgecolor='black')
    plt.title(f"Distribución de Distancias: {nombre}")
    plt.xlabel("Distancia (Å)"); plt.ylabel("Frecuencia")
    plt.savefig(os.path.join(ruta_carpeta, f"{nombre}_histograma.png"))
    plt.close()

def validar_estructuras(xyz_block):
    """Doble validación con OpenBabel y RDKit."""
    ob_ok, smi_ob = False, None
    rd_ok, smi_rd = False, None
    
    # 1. OpenBabel
    obConv = openbabel.OBConversion()
    obConv.SetInAndOutFormats("xyz", "smi")
    mol_ob = openbabel.OBMol()
    if obConv.ReadString(mol_ob, xyz_block):
        mol_ob.ConnectTheDots(); mol_ob.PerceiveBondOrders()
        smi_ob = obConv.WriteString(mol_ob).strip().split()[0]
        ob_ok = True

    # 2. RDKit
    try:
        mol_rd = Chem.MolFromXYZBlock(xyz_block)
        if mol_rd:
            rdDetermineBonds.DetermineBonds(mol_rd)
            smi_rd = Chem.MolToSmiles(mol_rd)
            rd_ok = True
    except: pass

    return ob_ok, rd_ok, smi_ob, smi_rd

def ejecutar_procesamiento():
    archivos = glob.glob(os.path.join(INPUT_DIR, "*.sumviz"))
    
    if not archivos:
        print(f"\n[!] No se encontraron archivos en: {INPUT_DIR}")
        return

    for ruta in archivos:
        # 1. Definir el nombre basado en el archivo .sumviz
        nombre = os.path.splitext(os.path.basename(ruta))[0]
        
        # 2. Crear la carpeta específica para esta molécula
        ruta_salida = os.path.join(OUTPUT_BASE_DIR, nombre)
        ruta_debug = os.path.join(ruta_salida, "debug")
        os.makedirs(ruta_debug, exist_ok=True)
        
        generar_readme(ruta_salida)
        reporte_auditoria = [] # Reporte individual por carpeta

        coords, puntos, atomos = extraer_datos(ruta)
        
        if coords:
            xyz_block = f"{len(coords)}\n{nombre}\n" + "\n".join(coords)
            matriz, lista_dist, dist_max, num_anomalias = calcular_matriz_y_max(puntos, atomos)
            
            ob_ok, rd_ok, smi_ob, smi_rd = validar_estructuras(xyz_block)
            
            es_mala = not (ob_ok and rd_ok) or (dist_max > 15.0) or (num_anomalias > 0)
            
            # El destino ahora es dentro de la carpeta con el nombre de la molécula
            destino = ruta_debug if es_mala else ruta_salida
            
            with open(os.path.join(destino, f"{nombre}.xyz"), "w") as f: f.write(xyz_block)
            with open(os.path.join(destino, f"{nombre}_distancias.json"), "w") as f: json.dump(matriz, f, indent=2)
            generar_histograma(lista_dist, nombre, destino)

            reporte_auditoria.append({
                "archivo": nombre,
                "status": "DEBUG" if es_mala else "PASSED",
                "distancia_max": round(dist_max, 4),
                "anomalias_vecinos": num_anomalias,
                "formula": dict(Counter(atomos)),
                "smiles_ob": smi_ob,
                "smiles_rd": smi_rd
            })
            
            # Guardar el JSON dentro de su carpeta correspondiente
            with open(os.path.join(ruta_salida, "auditoria_final.json"), "w") as f:
                json.dump(reporte_auditoria, f, indent=4)
                
            print(f"{'DEBUG' if es_mala else 'OK'}: Carpeta creada -> {nombre} ({num_anomalias})")

    print(f"\nProceso finalizado.")

if __name__ == "__main__":
    ejecutar_procesamiento()