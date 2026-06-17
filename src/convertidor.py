import os
import glob
import json
import math
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from openbabel import openbabel
from datetime import datetime
from collections import Counter
import tkinter as tk
from tkinter import messagebox

# --- CONFIGURACIÓN DE RUTAS BASE ---
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Generamos un timestamp del día actual para agrupar el experimento
FECHA_EXPERIMENTO = datetime.now().strftime("%Y-%m-%d")

# Directorio Raíz de Salidas organizado por Fecha
EXP_OUTPUT_DIR = os.path.join(BASE_DIR, "data", "experimentos", FECHA_EXPERIMENTO)

# Subcarpetas estándar para archivos aprobados (PASSED)
XYZ_DIR = os.path.join(EXP_OUTPUT_DIR, "coordenadas_xyz")
JSON_DIR = os.path.join(EXP_OUTPUT_DIR, "reportes_json")
PLOTS_DIR = os.path.join(EXP_OUTPUT_DIR, "plots_histogramas")
PLOTS_INDIV_DIR = os.path.join(PLOTS_DIR, "plots_individuales")

# --- NUEVO SISTEMA DE DEBUG INTEGRADO ---
# Subcarpetas físicas para el aislamiento de estructuras con errores o anomalías
DEBUG_BASE_DIR = os.path.join(EXP_OUTPUT_DIR, "debug")
DEBUG_XYZ_DIR = os.path.join(DEBUG_BASE_DIR, "coordenadas_xyz")
DEBUG_JSON_DIR = os.path.join(DEBUG_BASE_DIR, "reportes_json")
DEBUG_PLOTS_DIR = os.path.join(DEBUG_BASE_DIR, "plots_individuales")

# --- REFERENCIAS QUÍMICAS ---
RADIOS_COVALENTES = {
    "H": 0.32,  "He": 0.46,
    "Li": 1.33, "Be": 1.02, "B": 0.85,  "C": 0.75,  "N": 0.71,  "O": 0.63,  "F": 0.64,  "Ne": 0.67,
    "Na": 1.55, "Mg": 1.39, "Al": 1.26, "Si": 1.16, "P": 1.11, "S": 1.03, "Cl": 0.99, "Ar": 0.96,
    "K": 1.96,  "Ca": 1.71, "Br": 1.14, "I": 1.33
}

def preguntar_generar_histograma():
    """Muestra una ventana flotante para decidir si se procesan los gráficos."""
    root = tk.Tk()
    root.withdraw()
    root.attributes("-topmost", True)
    respuesta = messagebox.askyesno(
        title="Configuración de Gráficos",
        message=f"¿Deseas generar los histogramas (Individuales y Global) para el experimento de hoy ({FECHA_EXPERIMENTO})?"
    )
    root.destroy()
    return respuesta

def generar_readme(ruta_carpeta):
    """Crea la documentación técnica del proceso en la raíz del experimento."""
    contenido = f"""SISTEMA DE PROCESAMIENTO MOLECULAR + AUDITORÍA AVANZADA
==================================================
Fecha del Experimento : {FECHA_EXPERIMENTO}
Hora de Ejecución     : {datetime.now().strftime("%H:%M:%S")}
Estadísticas de Lote  : Organizado jerárquicamente en XYZ, JSON y PLOTS (Con aislamiento DEBUG).
--------------------------------------------------
- Filtro de Calidad: Validación por radios covalentes (choque dinámico adaptativo).
- Umbral Inferior Absoluto: Alerta si d < 0.55 Å (Límite de colisión física nuclear).
- Validación Cruzada: Comparativa SMILES OpenBabel vs RDKit.
- Test Visual: Histogramas individuales y globales de distribución de distancias con KDE adaptativo.
- Límites de Distancia: Alerta en distancias excesivas > 15.0 Å.
--------------------------------------------------
"""
    with open(os.path.join(ruta_carpeta, "README_PROCESO.txt"), "w", encoding="utf-8") as f:
        f.write(contenido)

def extraer_datos(ruta_archivo):
    """Extrae coordenadas y convierte de Bohr a Ángstrom manejando encabezados de texto."""
    coordenadas, puntos, atomos = [], [], []
    factor = 0.529177 
    en_seccion = False

    try:
        with open(ruta_archivo, 'r', encoding='utf-8') as f:
            for linea in f:
                if "Nuclear Charges and Cartesian Coordinates" in linea:
                    en_seccion = True; continue
                if en_seccion:
                    partes = linea.split()
                    if len(partes) == 5:
                        simbolo = "".join([c for c in partes[0] if c.isalpha()])
                        try:
                            x = float(partes[2]) * factor
                            y = float(partes[3]) * factor
                            z = float(partes[4]) * factor
                            
                            coordenadas.append(f"{simbolo} {x:10.6f} {y:10.6f} {z:10.6f}")
                            puntos.append((x, y, z))
                            atomos.append(simbolo)
                        except ValueError:
                            continue
                    elif len(coordenadas) > 0: 
                        break
        return coordenadas, puntos, atomos
    except Exception as e:
        print(f"Error leyendo {ruta_archivo}: {e}")
        return None, None, None

def auditar_distancias(puntos, atomos):
    """Calcula distancias, detecta colisiones reales y valida límites físicos."""
    matriz = []
    lista_dist = []
    anomalias = []
    n = len(puntos)
    
    for i in range(n):
        for j in range(i + 1, n):
            d = math.sqrt(sum((a - b) ** 2 for a, b in zip(puntos[i], puntos[j])))
            
            r1 = RADIOS_COVALENTES.get(atomos[i], 0.75)
            r2 = RADIOS_COVALENTES.get(atomos[j], 0.75)
            umbral_choque = (r1 + r2) * 0.55
            
            if d < 0.55 or d < umbral_choque:
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

def generar_histograma_individual(distancias, serial_id, ruta_destino):
    """Genera un análisis estadístico visual específico para un archivo del lote."""
    if not distancias: return
    
    min_real = min(distancias)
    max_real = max(distancias)
    
    plt.figure(figsize=(10, 5))
    
    if max_real < 15.0:
        plt.xlim(0, max_real + 1.0)
    else:
        plt.xlim(0, max_real + 2.0)
        
    sns.histplot(distancias, bins=25, kde=True, color='royalblue', edgecolor='black', alpha=0.6)
    
    plt.axvline(0.55, color='red', linestyle='--', linewidth=1.5, label='Límite Choque Inferior (0.55Å)')
    if max_real >= 14.0:
        plt.axvline(15.0, color='darkred', linestyle=':', linewidth=2, label='Límite Superior Máx (15.0Å)')
        
    plt.axvline(min_real, color='blue', linestyle='-.', linewidth=1.0, label=f'Mínimo Real ({min_real:.3f}Å)')
    plt.axvline(max_real, color='orange', linestyle='-.', linewidth=1.0, label=f'Máximo Real ({max_real:.3f}Å)')
    
    plt.title(f"Auditoría Individual de Distancias: {serial_id}")
    plt.xlabel("Distancia (Å)")
    plt.ylabel("Frecuencia")
    plt.legend(loc='upper right')
    plt.grid(axis='y', alpha=0.3)
    
    plt.savefig(os.path.join(ruta_destino, f"plot_{serial_id}.png"), bbox_inches='tight', dpi=150)
    plt.close()

def generar_histograma_global(distancias_globales, ruta_destino):
    """Genera un análisis estadístico visual enfocado con límites adaptativos inteligentes."""
    if not distancias_globales: 
        print("[!] No hay datos de distancias para graficar.")
        return
        
    min_real = min(distancias_globales)
    max_real = max(distancias_globales)
    
    plt.figure(figsize=(12, 6))
    
    if max_real < 15.0:
        plt.xlim(0, max_real + 1.0)
    else:
        plt.xlim(0, max_real + 2.0)
    
    sns.histplot(distancias_globales, bins=60, kde=True, color='purple', edgecolor='black', alpha=0.7)
    
    plt.axvline(0.55, color='red', linestyle='--', linewidth=1.5, label='Umbral Choque Inferior (0.55Å)')
    if max_real >= 14.0: 
        plt.axvline(15.0, color='darkred', linestyle=':', linewidth=2, label='Umbral Max Superior (15.0Å)')
    
    plt.axvline(min_real, color='blue', linestyle='-.', linewidth=1.2, label=f'Mínimo Real ({min_real:.3f}Å)')
    plt.axvline(max_real, color='orange', linestyle='-.', linewidth=1.2, label=f'Máximo Real ({max_real:.3f}Å)')
    
    hora_actual = datetime.now().strftime("%H:%M:%S")
    plt.title(f"Auditoría Global de Distancias - Experimento: {FECHA_EXPERIMENTO} ({hora_actual})", fontsize=13, fontweight='bold')
    plt.xlabel("Distancia Interatómica (Å)", fontsize=11)
    plt.ylabel("Frecuencia Acumulada", fontsize=11)
    plt.legend(loc='upper right', frameon=True, facecolor='white', edgecolor='none')
    plt.grid(axis='y', alpha=0.3)
    
    stamp_archivo = datetime.now().strftime("%H%M%S")
    nombre_archivo = f"plot_global_{stamp_archivo}.png"
    archivo_salida = os.path.join(ruta_destino, nombre_archivo)
    
    plt.savefig(archivo_salida, bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"\n[+] Histograma GLOBAL guardado en PLOTS -> {archivo_salida}")

def validar_estructuras(xyz_block):
    """Doble validación de conectividad química."""
    ob_ok, smi_ob = False, None
    rd_ok, smi_rd = False, None
    
    obConv = openbabel.OBConversion()
    obConv.SetInAndOutFormats("xyz", "smi")
    mol_ob = openbabel.OBMol()
    if obConv.ReadString(mol_ob, xyz_block):
        mol_ob.ConnectTheDots()
        mol_ob.PerceiveBondOrders()
        smi_ob = obConv.WriteString(mol_ob).strip().split()[0]
        ob_ok = True

    try:
        mol_rd = Chem.MolFromXYZBlock(xyz_block)
        if mol_rd:
            rdDetermineBonds.DetermineBonds(mol_rd)
            smi_rd = Chem.MolToSmiles(mol_rd)
            rd_ok = True
    except: pass

    return ob_ok, rd_ok, smi_ob, smi_rd

def ejecutar_procesamiento(rutas_input, activar_histograma):
    """Procesa los archivos .sumviz mapeándolos y separándolos físicamente por estatus."""
    archivos = []
    distancias_totales_globales = []
    
    for ruta_item in rutas_input:
        ruta_absolute = os.path.normpath(os.path.join(BASE_DIR, ruta_item)) if not os.path.isabs(ruta_item) else ruta_item

        if os.path.isfile(ruta_absolute):
            folder_objetivo = os.path.dirname(ruta_absolute)
            archivos_encontrados = glob.glob(os.path.join(folder_objetivo, "*.sumviz"))
            archivos.extend(archivos_encontrados)
        elif os.path.isdir(ruta_absolute):
            archivos_encontrados = glob.glob(os.path.join(ruta_absolute, "*.sumviz"))
            archivos.extend(archivos_encontrados)
            print(f"[+] Buscando en carpeta: {ruta_absolute} ({len(archivos_encontrados)} archivos)")

    archivos = list(dict.fromkeys(archivos))
    if not archivos:
        print(f"\n[!] Sin archivos .sumviz para procesar.")
        return

    # Crear directorios estándar
    os.makedirs(XYZ_DIR, exist_ok=True)
    os.makedirs(JSON_DIR, exist_ok=True)
    os.makedirs(PLOTS_DIR, exist_ok=True)
    if activar_histograma:
        os.makedirs(PLOTS_INDIV_DIR, exist_ok=True)
        
    generar_readme(EXP_OUTPUT_DIR)

    print(f"\n[*] Procesando lote en el experimento {FECHA_EXPERIMENTO} ({len(archivos)} archivos)...")

    for idx, ruta in enumerate(archivos, start=1):
        nombre = os.path.splitext(os.path.basename(ruta))[0]
        serial_id = f"{idx:03d}_{nombre}"
        
        coords, puntos, atomos = extraer_datos(ruta)
        
        if coords:
            xyz_block = f"{len(coords)}\n{nombre}\n" + "\n".join(coords)
            matriz, lista_dist, dist_max, anomalias = auditar_distancias(puntos, atomos)
            ob_ok, rd_ok, smi_ob, smi_rd = validar_estructuras(xyz_block)
            
            distancias_totales_globales.extend(lista_dist)
            num_choques = len(anomalias)
            
            # Evaluación lógica de calidad
            es_anomalo = not (ob_ok and rd_ok) or (dist_max > 15.0) or (num_choques > 0)
            status_str = "DEBUG" if es_anomalo else "PASSED"
            
            # REDIRECCIÓN DINÁMICA DE RUTAS: Si falla, crea y asigna las carpetas del bloque debug
            if es_anomalo:
                os.makedirs(DEBUG_XYZ_DIR, exist_ok=True)
                os.makedirs(DEBUG_JSON_DIR, exist_ok=True)
                if activar_histograma: os.makedirs(DEBUG_PLOTS_DIR, exist_ok=True)
                
                destino_xyz = DEBUG_XYZ_DIR
                destino_json = DEBUG_JSON_DIR
                destino_plots = DEBUG_PLOTS_DIR
            else:
                destino_xyz = XYZ_DIR
                destino_json = JSON_DIR
                destino_plots = PLOTS_INDIV_DIR

            # --- GUARDADO EN DIRECTORIO ASIGNADO ---
            # 1. Guardar Estructura Molecular (.xyz)
            with open(os.path.join(destino_xyz, f"{serial_id}.xyz"), "w") as f: 
                f.write(xyz_block)
            
            # 2. Guardar Reportes Estructurados (.json)
            reporte_distancias = {"matriz": matriz, "alertas": anomalias}
            with open(os.path.join(destino_json, f"{serial_id}_distancias.json"), "w") as f: 
                json.dump(reporte_distancias, f, indent=2)

            auditoria_molecular = {
                "identificador_serial": serial_id,
                "archivo_origen": nombre,
                "fecha_procesamiento": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "status_calidad": status_str,
                "metrics": {
                    "dist_max_A": round(dist_max, 4),
                    "num_anomalias": num_choques,
                    "formula_estequiometrica": dict(Counter(atomos))
                },
                "quimica": {
                    "smiles_openbabel": smi_ob,
                    "smiles_rdkit": smi_rd,
                    "validacion_cruzada_smiles": (smi_ob == smi_rd)
                },
                "detalles_anomalias": anomalias
            }
            
            with open(os.path.join(destino_json, f"{serial_id}_auditoria.json"), "w") as f:
                json.dump(auditoria_molecular, f, indent=4)
            
            # 3. Guardar el histograma individual en el destino correspondiente
            if activar_histograma:
                generar_histograma_individual(lista_dist, serial_id, destino_plots)
                
            print(f" [{status_str}]: {serial_id} | Choques: {num_choques} | Max: {dist_max:.2f}Å")

    # --- GENERACIÓN DEL PLOT GLOBAL ---
    if activar_histograma:
        generar_histograma_global(distancias_totales_globales, PLOTS_DIR)

    print(f"\n--- Experimento {FECHA_EXPERIMENTO} Finalizado con Éxito ---")
    print(f"[➔] Datos estructurados en: {EXP_OUTPUT_DIR}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Procesador Molecular con Serialización Cronológica Estricta.")
    parser.add_argument('rutas', nargs='*', help='Rutas de carpetas con archivos .sumviz')
    args = parser.parse_args()
    
    hacer_histograma = preguntar_generar_histograma()
    
    if not args.rutas:
        ruta_por_defecto = os.path.join(BASE_DIR, "data", "input")
        print(f"[*] Usando ruta por defecto: {ruta_por_defecto}")
        ejecutar_procesamiento([ruta_por_defecto], hacer_histograma)
    else:
        ejecutar_procesamiento(args.rutas, hacer_histograma)