import os
import glob
import json
import math
import shutil
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from openbabel import openbabel
from datetime import datetime
from collections import Counter

# --- CONFIGURACIÓN DE RUTAS ---
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUT_DIR = os.path.join(BASE_DIR, "data", "input")
OUTPUT_BASE_DIR = os.path.join(BASE_DIR, "data", "output")

# --- REFERENCIAS QUÍMICAS ---
# Radios covalentes aproximados en Ångstroms para validación de proximidad
RADIOS_COVALENTES = {
    "H": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96, "B": 0.84,
    "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58,
    "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07,
    "S": 1.05, "Cl": 1.02, "Ar": 1.06, "K": 2.03, "Ca": 1.76
}

def generar_readme(ruta_carpeta):
    """Crea la documentación técnica del proceso."""
    contenido = f"""
SISTEMA DE PROCESAMIENTO MOLECULAR + AUDITORÍA AVANZADA
Fecha: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
--------------------------------------------------
- Filtro de Calidad: Validación por radios covalentes (choque atómico).
- Validación Cruzada: Comparativa SMILES OpenBabel vs RDKit.
- Test Visual: Histograma de distribución de distancias con KDE.
- Límites de Distancia: Alerta en distancias < 0.8 Å o > 15.0 Å.
--------------------------------------------------
"""
    with open(os.path.join(ruta_carpeta, "README_PROCESO.txt"), "w", encoding="utf-8") as f:
        f.write(contenido)

def extraer_datos(ruta_archivo):
    """Extrae coordenadas y convierte de Bohr a Ángstrom."""
    coordenadas, puntos, atomos = [], [], []
    factor = 0.529177 # Bohr to Angstrom
    en_seccion = False
    elementos_validos = set(RADIOS_COVALENTES.keys())

    try:
        with open(ruta_archivo, 'r', encoding='utf-8') as f:
            for linea in f:
                if "Nuclear Charges and Cartesian Coordinates" in linea:
                    en_seccion = True; continue
                if en_seccion:
                    partes = linea.split()
                    if len(partes) == 5:
                        simbolo = "".join([c for c in partes[0] if c.isalpha()])
                        if simbolo in elementos_validos:
                            x, y, z = float(partes[2])*factor, float(partes[3])*factor, float(partes[4])*factor
                            coordenadas.append(f"{simbolo} {x:10.6f} {y:10.6f} {z:10.6f}")
                            puntos.append((x, y, z))
                            atomos.append(simbolo)
                    elif len(coordenadas) > 0: break
        return coordenadas, puntos, atomos
    except Exception as e:
        print(f"Error leyendo {ruta_archivo}: {e}")
        return None, None, None

def auditar_distancias(puntos, atomos):
    """Calcula distancias, detecta colisiones y valida límites físicos."""
    matriz = []
    lista_dist = []
    anomalias = []
    n = len(puntos)
    
    for i in range(n):
        for j in range(i + 1, n):
            d = math.sqrt(sum((a - b) ** 2 for a, b in zip(puntos[i], puntos[j])))
            
            # Validación de choque (basado en radios covalentes)
            r1 = RADIOS_COVALENTES.get(atomos[i], 0.7)
            r2 = RADIOS_COVALENTES.get(atomos[j], 0.7)
            umbral_choque = (r1 + r2) * 0.6 # 60% de la suma de radios es el límite crítico
            
            if d < 0.8 or d < umbral_choque:
                anomalias.append({
                    "par": f"{atomos[i]}({i})-{atomos[j]}({j})",
                    "distancia": round(d, 4),
                    "tipo": "Choque Atómico"
                })
            
            matriz.append({
                "de": f"{atomos[i]}({i})", 
                "a": f"{atomos[j]}({j})", 
                "dist_A": round(d, 4)
            })
            lista_dist.append(d)
            
    dist_max = max(lista_dist) if lista_dist else 0
    return matriz, lista_dist, dist_max, anomalias

def generar_histograma(distancias, nombre, ruta_carpeta):
    """Genera análisis estadístico visual de las distancias."""
    if not distancias: return
    plt.figure(figsize=(9, 5))
    
    # Usar seaborn para una distribución más estética
    sns.histplot(distancias, bins=15, kde=True, color='royalblue', edgecolor='black')
    
    # Línea de advertencia para distancias cortas
    plt.axvline(0.8, color='red', linestyle='--', label='Límite Choque (0.8Å)')
    
    plt.title(f"Distribución de Distancias: {nombre}")
    plt.xlabel("Distancia (Å)")
    plt.ylabel("Frecuencia")
    plt.legend()
    plt.grid(axis='y', alpha=0.3)
    
    plt.savefig(os.path.join(ruta_carpeta, f"{nombre}_distribucion.png"))
    plt.close()

def validar_estructuras(xyz_block):
    """Doble validación de conectividad química."""
    ob_ok, smi_ob = False, None
    rd_ok, smi_rd = False, None
    
    # OpenBabel
    obConv = openbabel.OBConversion()
    obConv.SetInAndOutFormats("xyz", "smi")
    mol_ob = openbabel.OBMol()
    if obConv.ReadString(mol_ob, xyz_block):
        mol_ob.ConnectTheDots()
        mol_ob.PerceiveBondOrders()
        smi_ob = obConv.WriteString(mol_ob).strip().split()[0]
        ob_ok = True

    # RDKit
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
        print(f"\n[!] No hay archivos .sumviz en: {INPUT_DIR}")
        return

    for ruta in archivos:
        nombre = os.path.splitext(os.path.basename(ruta))[0]
        ruta_salida = os.path.join(OUTPUT_BASE_DIR, nombre)
        ruta_debug = os.path.join(ruta_salida, "debug")
        os.makedirs(ruta_debug, exist_ok=True)
        
        generar_readme(ruta_salida)
        coords, puntos, atomos = extraer_datos(ruta)
        
        if coords:
            xyz_block = f"{len(coords)}\n{nombre}\n" + "\n".join(coords)
            matriz, lista_dist, dist_max, anomalias = auditar_distancias(puntos, atomos)
            ob_ok, rd_ok, smi_ob, smi_rd = validar_estructuras(xyz_block)
            
            # CRITERIO DE CALIDAD:
            # Pasa si: Ambos validan SMILES Y dist_max < 15Å Y no hay choques atómicos
            num_choques = len(anomalias)
            es_mala = not (ob_ok and rd_ok) or (dist_max > 15.0) or (num_choques > 0)
            
            destino = ruta_debug if es_mala else ruta_salida
            
            # Guardar resultados
            with open(os.path.join(destino, f"{nombre}.xyz"), "w") as f: f.write(xyz_block)
            with open(os.path.join(destino, f"{nombre}_distancias.json"), "w") as f: 
                json.dump({"matriz": matriz, "alertas": anomalias}, f, indent=2)
            
            generar_histograma(lista_dist, nombre, destino)

            # Reporte final detallado
            auditoria = {
                "archivo": nombre,
                "status": "DEBUG" if es_mala else "PASSED",
                "metrics": {
                    "dist_max": round(dist_max, 4),
                    "num_anomalias": num_choques,
                    "formula": dict(Counter(atomos))
                },
                "quimica": {
                    "smiles_ob": smi_ob,
                    "smiles_rd": smi_rd,
                    "valida_cruzada": (smi_ob == smi_rd)
                },
                "lista_anomalias": anomalias
            }
            
            with open(os.path.join(ruta_salida, "auditoria_final.json"), "w") as f:
                json.dump(auditoria, f, indent=4)
                
            status_str = "DEBUG" if es_mala else "PASSED"
            print(f"{status_str}: {nombre} | Choques: {num_choques} | Max: {dist_max:.2f}Å")

    print(f"\n--- Proceso finalizado ---")

if __name__ == "__main__":
    ejecutar_procesamiento()