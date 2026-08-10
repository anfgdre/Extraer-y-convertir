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

# Referencia orientativa de valencia/vecinos esperados para determinar si un átomo está "rodeado/completo"
VALENCIA_REFERENCIA = {
    "H": 1,  "F": 1,  "Cl": 1, "Br": 1, "I": 1,
    "O": 2,  "S": 2,
    "N": 3,  "P": 3,  "B": 3,
    "C": 4,  "Si": 4, "Ge": 4
}

def preguntar_generar_reportes_y_graficos():
    """Muestra una ventana flotante para decidir si se procesan los gráficos y reportes individuales exitosos."""
    root = tk.Tk()
    root.withdraw()
    root.attributes("-topmost", True)
    respuesta = messagebox.askyesno(
        title="Configuración de Reportes y Gráficos",
        message=f"¿Deseas generar los reportes individuales JSON e histogramas para los casos exitosos (PASSED) del experimento de hoy ({FECHA_EXPERIMENTO})?\n\n(Nota: Los casos con anomalías se generarán siempre en DEBUG)"
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
- Filtro de Calidad: Auditoría de distancias mínimas por par (vecino más cercano).
- Umbrales de Choque: Alerta si d < 0.95 Å y Choque Crítico si d < 0.80 Å (átomos encimados).
- Validación Cruzada Triple: Matriz de distancias por pares geométrica vs OpenBabel vs RDKit.
- Test Visual: Histogramas enfocados en distancias cortas (0.0 Å a 2.5 Å) con KDE adaptativo.
- Análisis de Entorno: Conexiones reales por átomo, detección de fragmentos huérfanos y estado de coordinación.
--------------------------------------------------
"""
    os.makedirs(ruta_carpeta, exist_ok=True)
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

def auditar_distancias_minimas(puntos, atomos, umbral_alerta_encimado=0.95, umbral_choque_critico=0.80):
    """
    Obtiene la matriz de distancias por pares, recupera la distancia MÍNIMA 
    de cada átomo a su vecino más cercano y detecta átomos encimados.
    
    Retorna:
    - distancias_minimas_por_atomo: Lista con la menor distancia registrada para cada átomo.
    - distancias_cortas_par: Distancias en el rango de interacción cercana (0.5 Å - 2.5 Å).
    - anomalias: Átomos detectados con solapamiento/choque o distancia anormalmente corta.
    - min_global: La distancia mínima absoluta en toda la estructura.
    """
    n = len(puntos)
    if n < 2:
        return [], [], [], 0.0, []

    matriz_dist = [[0.0] * n for _ in range(n)]
    matriz_reporte = []
    distancias_minimas_por_atomo = []
    distancias_cortas_par = []
    anomalias = []
    pares_evaluados = set()

    # 1. Matriz completa de distancias por pares
    for i in range(n):
        for j in range(i + 1, n):
            d = math.sqrt(sum((a - b) ** 2 for a, b in zip(puntos[i], puntos[j])))
            matriz_dist[i][j] = d
            matriz_dist[j][i] = d

            matriz_reporte.append({
                "de": f"{atomos[i]}({i})", 
                "a": f"{atomos[j]}({j})", 
                "dist_A": round(d, 4)
            })

            if 0.5 <= d <= 2.5:
                distancias_cortas_par.append(round(d, 4))

    # 2. Búsqueda del vecino más cercano para cada átomo
    for i in range(n):
        distancias_vecinos = [(matriz_dist[i][j], j) for j in range(n) if i != j]
        
        if not distancias_vecinos:
            continue
            
        dist_min, idx_vecino = min(distancias_vecinos, key=lambda x: x[0])
        distancias_minimas_por_atomo.append(round(dist_min, 4))

        r1 = RADIOS_COVALENTES.get(atomos[i], 0.75)
        r2 = RADIOS_COVALENTES.get(atomos[idx_vecino], 0.75)
        umbral_dinamico = (r1 + r2) * 0.60  # Factor de choque físico relativo

        # Evaluamos anomalías sin duplicar pares
        par_key = tuple(sorted([i, idx_vecino]))
        if par_key not in pares_evaluados:
            if dist_min < umbral_choque_critico or dist_min < umbral_dinamico:
                anomalias.append({
                    "par": f"{atomos[i]}({i})-{atomos[idx_vecino]}({idx_vecino})",
                    "distancia_A": round(dist_min, 4),
                    "tipo_error": "Choque Crítico / Átomos Encimados",
                    "umbral_corte_A": round(umbral_dinamico, 4)
                })
                pares_evaluados.add(par_key)
            elif dist_min < umbral_alerta_encimado:
                anomalias.append({
                    "par": f"{atomos[i]}({i})-{atomos[idx_vecino]}({idx_vecino})",
                    "distancia_A": round(dist_min, 4),
                    "tipo_error": f"Alerta Distancia Corta (< {umbral_alerta_encimado} Å)",
                    "umbral_corte_A": umbral_alerta_encimado
                })
                pares_evaluados.add(par_key)

    min_global = min(distancias_minimas_por_atomo) if distancias_minimas_por_atomo else 0.0

    return distancias_minimas_por_atomo, distancias_cortas_par, anomalias, min_global, matriz_reporte

def generar_histograma_distancias_cortas(distancias_minimas, serial_id, ruta_destino):
    """
    Grafica el histograma de las distancias MÍNIMAS al vecino más cercano,
    enfocado estrictamente en el rango relevante para control de calidad (0.0 Å - 2.5 Å).
    """
    if not distancias_minimas: 
        return
    
    min_real = min(distancias_minimas)
    max_real = max(distancias_minimas)
    
    plt.figure(figsize=(9, 4.5))
    
    plt.xlim(0.0, max(2.5, max_real + 0.2))
        
    sns.histplot(
        distancias_minimas, 
        bins=25, 
        kde=True, 
        color='crimson' if min_real < 0.80 else 'royalblue', 
        edgecolor='black', 
        alpha=0.6
    )
    
    plt.axvline(0.80, color='red', linestyle='--', linewidth=1.5, label='Choque Crítico (0.80Å)')
    plt.axvline(0.95, color='orange', linestyle=':', linewidth=1.5, label='Alerta Encimado (0.95Å)')
    plt.axvline(min_real, color='darkblue', linestyle='-.', linewidth=1.2, label=f'Mínimo Real ({min_real:.3f}Å)')
    
    plt.title(f"Auditoría de Distancias Mínimas al Vecino Más Cercano: {serial_id}", fontsize=11, fontweight='bold')
    plt.xlabel("Distancia al Vecino Más Cercano (Å)", fontsize=10)
    plt.ylabel("Frecuencia de Átomos", fontsize=10)
    plt.legend(loc='upper right', fontsize=8)
    plt.grid(axis='y', alpha=0.3)
    
    os.makedirs(ruta_destino, exist_ok=True)
    plt.savefig(os.path.join(ruta_destino, f"plot_{serial_id}.png"), bbox_inches='tight', dpi=150)
    plt.close()

def obtener_vecinos_por_radio(puntos, atomos, idx_atomo=1, radio_angstrom=2.0):
    """Retorna los átomos dentro de un radio fijo en Ångstroms alrededor de un átomo objetivo."""
    p_ref = puntos[idx_atomo]
    vecinos_en_radio = []
    
    for j, p_o in enumerate(puntos):
        if j == idx_atomo:
            continue
        d = math.sqrt(sum((a - b) ** 2 for a, b in zip(p_ref, p_o)))
        if d <= radio_angstrom:
            vecinos_en_radio.append({
                "idx": j,
                "simbolo": atomos[j],
                "distancia_A": round(d, 4)
            })
            
    return vecinos_en_radio, len(vecinos_en_radio)

def generar_histograma_global(distancias_minimas_globales, ruta_destino):
    """Genera un análisis estadístico acumulado del lote enfocado en distancias cortas/mínimas."""
    if not distancias_minimas_globales: 
        print("[!] No hay datos de distancias para graficar.")
        return
        
    min_real = min(distancias_minimas_globales)
    max_real = max(distancias_minimas_globales)
    
    plt.figure(figsize=(11, 5.5))
    plt.xlim(0.0, max(2.5, max_real + 0.2))
    
    sns.histplot(
        distancias_minimas_globales, 
        bins=50, 
        kde=True, 
        color='purple', 
        edgecolor='black', 
        alpha=0.7
    )
    
    plt.axvline(0.80, color='red', linestyle='--', linewidth=1.5, label='Límite Encimado Crítico (0.80Å)')
    plt.axvline(0.95, color='orange', linestyle=':', linewidth=1.5, label='Umbral Alerta Encimado (0.95Å)')
    plt.axvline(min_real, color='blue', linestyle='-.', linewidth=1.2, label=f'Mínimo Global ({min_real:.3f}Å)')
    
    hora_actual = datetime.now().strftime("%H:%M:%S")
    plt.title(f"Auditoría Global de Distancias Mínimas - Experimento: {FECHA_EXPERIMENTO} ({hora_actual})", fontsize=12, fontweight='bold')
    plt.xlabel("Distancia al Vecino Más Cercano (Å)", fontsize=10)
    plt.ylabel("Frecuencia Acumulada", fontsize=10)
    plt.legend(loc='upper right', frameon=True, facecolor='white')
    plt.grid(axis='y', alpha=0.3)
    
    stamp_archivo = datetime.now().strftime("%H%M%S")
    nombre_archivo = f"plot_global_{stamp_archivo}.png"
    
    os.makedirs(ruta_destino, exist_ok=True)
    archivo_salida = os.path.join(ruta_destino, nombre_archivo)
    
    plt.savefig(archivo_salida, bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"\n[+] Histograma GLOBAL guardado en PLOTS -> {archivo_salida}")

def determinar_conectividad_por_pares(puntos, atomos, factor_tolerancia=1.20):
    """Determina conectividad por radios covalentes escalados por factor de tolerancia."""
    n = len(atomos)
    vecinos_por_pares = {i: [] for i in range(n)}
    
    for i in range(n):
        r1 = RADIOS_COVALENTES.get(atomos[i], 0.75)
        for j in range(i + 1, n):
            r2 = RADIOS_COVALENTES.get(atomos[j], 0.75)
            umbral_enlace = (r1 + r2) * factor_tolerancia
            
            d = math.sqrt(sum((a - b) ** 2 for a, b in zip(puntos[i], puntos[j])))
            
            if 0.55 <= d <= umbral_enlace:
                vecinos_por_pares[i].append(j)
                vecinos_por_pares[j].append(i)
                
    for i in vecinos_por_pares:
        vecinos_por_pares[i].sort()
        
    return vecinos_por_pares

def validar_y_analizar_entorno(xyz_block, puntos, atomos):
    """
    Auditoría estricta comparando:
    1. Distancias directas por pares (Geometría pura)
    2. Enlaces percibidos por RDKit
    3. Enlaces percibidos por OpenBabel
    """
    ob_ok, smi_ob = False, None
    rd_ok, smi_rd = False, None
    analisis_entorno = []
    tiene_atomos_aislados = False
    coinciden_vecinos_global = True
    
    vecinos_pares = determinar_conectividad_por_pares(puntos, atomos, factor_tolerancia=1.20)
    
    vecinos_openbabel = {}
    obConv = openbabel.OBConversion()
    obConv.SetInAndOutFormats("xyz", "smi")
    mol_ob = openbabel.OBMol()
    if obConv.ReadString(mol_ob, xyz_block):
        mol_ob.ConnectTheDots()
        mol_ob.PerceiveBondOrders()
        smi_ob = obConv.WriteString(mol_ob).strip().split()[0]
        ob_ok = True
        
        for atomo_ob in openbabel.OBMolAtomIter(mol_ob):
            idx_0 = atomo_ob.GetIdx() - 1
            vecinos_idx = [vecino_ob.GetIdx() - 1 for vecino_ob in openbabel.OBAtomAtomIter(atomo_ob)]
            vecinos_openbabel[idx_0] = sorted(vecinos_idx)

    try:
        mol_rd = Chem.MolFromXYZBlock(xyz_block)
        if mol_rd:
            rdDetermineBonds.DetermineBonds(mol_rd)
            smi_rd = Chem.MolToSmiles(mol_rd)
            rd_ok = True
            
            for atomo in mol_rd.GetAtoms():
                idx = atomo.GetIdx()
                simbolo = atomo.GetSymbol()
                
                vecinos_rd_obj = atomo.GetNeighbors()
                vecinos_rd_idx = sorted([v.GetIdx() for v in vecinos_rd_obj])
                vecinos_rd_nombres = [f"{v.GetSymbol()}({v.GetIdx()})" for v in vecinos_rd_obj]
                
                vecinos_ob_idx = vecinos_openbabel.get(idx, [])
                vecinos_par_idx = vecinos_pares.get(idx, [])
                
                num_vecinos = len(vecinos_par_idx)
                if num_vecinos == 0:
                    tiene_atomos_aislados = True
                
                coinciden_vecinos_atomo = (vecinos_rd_idx == vecinos_ob_idx == vecinos_par_idx)
                if not coinciden_vecinos_atomo:
                    coinciden_vecinos_global = False
                
                max_esperado = VALENCIA_REFERENCIA.get(simbolo, 4)
                esta_rodeado = num_vecinos >= max_esperado
                
                analisis_entorno.append({
                    "atomo_idx": idx,
                    "simbolo": simbolo,
                    "cantidad_vecinos_geometria_pares": num_vecinos,
                    "cantidad_vecinos_rdkit": len(vecinos_rd_idx),
                    "cantidad_vecinos_openbabel": len(vecinos_ob_idx),
                    "vecinos_pares_idx": vecinos_par_idx,
                    "vecinos_rd_idx": vecinos_rd_idx,
                    "vecinos_ob_idx": vecinos_ob_idx,
                    "coincidencia_exacta_metodos": coinciden_vecinos_atomo,
                    "vecinos_conectados_rdkit": vecinos_rd_nombres,
                    "esta_rodeado": esta_rodeado,
                    "status_coordinacion": "AISLADO / HUÉRFANO" if num_vecinos == 0 else ("RODEADO / COMPLETO" if esta_rodeado else "CONECTADO PARCIAL")
                })
    except Exception as e: 
        coinciden_vecinos_global = False

    return ob_ok, rd_ok, smi_ob, smi_rd, analisis_entorno, tiene_atomos_aislados, coinciden_vecinos_global

def ejecutar_procesamiento(rutas_input, activar_reportes_opcionales):
    """Procesa los archivos .sumviz mapeándolos y separándolos físicamente por estatus."""
    archivos = []
    distancias_minimas_globales = []
    
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

    generar_readme(EXP_OUTPUT_DIR)

    print(f"\n[*] Procesando lote en el experimento {FECHA_EXPERIMENTO} ({len(archivos)} archivos)...")

    for idx, ruta in enumerate(archivos, start=1):
        nombre = os.path.splitext(os.path.basename(ruta))[0]
        serial_id = f"{idx:03d}_{nombre}"
        
        coords, puntos, atomos = extraer_datos(ruta)
        
        if coords:
            xyz_block = f"{len(coords)}\n{nombre}\n" + "\n".join(coords)
            
            # Auditoría enfocada en MÍNIMAS DISTANCIAS y SOLAPAMIENTOS
            dist_minimas_atomos, dist_cortas_par, anomalias_choque, min_dist_global, matriz_reporte = auditar_distancias_minimas(
                puntos, atomos, umbral_alerta_encimado=0.95, umbral_choque_critico=0.80
            )
            
            distancias_minimas_globales.extend(dist_minimas_atomos)
            
            ob_ok, rd_ok, smi_ob, smi_rd, analisis_entorno, tiene_atomos_aislados, coinciden_vecinos_global = validar_y_analizar_entorno(xyz_block, puntos, atomos)
            
            num_choques = len(anomalias_choque)
            
            # Evaluación lógica de calidad
            es_anomalo = (min_dist_global < 0.80) or (num_choques > 0) or not (ob_ok and rd_ok) or tiene_atomos_aislados or not coinciden_vecinos_global
            status_str = "DEBUG" if es_anomalo else "PASSED"
            
            debe_escribir_reporte = es_anomalo or activar_reportes_opcionales
            
            if es_anomalo:
                destino_xyz = DEBUG_XYZ_DIR
                destino_json = DEBUG_JSON_DIR
                destino_plots = DEBUG_PLOTS_DIR
            else:
                destino_xyz = XYZ_DIR
                destino_json = JSON_DIR
                destino_plots = PLOTS_INDIV_DIR

            # 1. Guardar Estructura Molecular (.xyz)
            os.makedirs(destino_xyz, exist_ok=True)
            with open(os.path.join(destino_xyz, f"{serial_id}.xyz"), "w") as f: 
                f.write(xyz_block)
            
            # 2. Guardar Reportes Estructurados (.json)
            if debe_escribir_reporte:
                os.makedirs(destino_json, exist_ok=True)
                
                reporte_distancias = {
                    "distancias_minimas_por_atomo": dist_minimas_atomos,
                    "distancia_minima_absoluta_A": min_dist_global,
                    "matriz_por_pares": matriz_reporte,
                    "alertas_solapamiento": anomalias_choque
                }
                with open(os.path.join(destino_json, f"{serial_id}_distancias.json"), "w") as f: 
                    json.dump(reporte_distancias, f, indent=2)

                auditoria_molecular = {
                    "identificador_serial": serial_id,
                    "archivo_origen": nombre,
                    "fecha_procesamiento": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "status_calidad": status_str,
                    "metrics": {
                        "dist_minima_absoluta_A": round(min_dist_global, 4),
                        "num_anomalias_choque": num_choques,
                        "tiene_atomos_aislados": tiene_atomos_aislados,
                        "formula_estequiometrica": dict(Counter(atomos))
                    },
                    "quimica": {
                        "smiles_openbabel": smi_ob,
                        "smiles_rdkit": smi_rd,
                        "validacion_cruzada_smiles": (smi_ob == smi_rd),
                        "coincidencia_exacta_metodos_conectividad": coinciden_vecinos_global
                    },
                    "auditoria_entornos_conectados_geometria_vs_rdkit_vs_openbabel": analisis_entorno,
                    "detalles_anomalias": anomalias_choque
                }
                
                with open(os.path.join(destino_json, f"{serial_id}_auditoria.json"), "w") as f:
                    json.dump(auditoria_molecular, f, indent=4)
            
            # 3. Guardar el histograma individual de distancias mínimas
            if debe_escribir_reporte:
                generar_histograma_distancias_cortas(dist_minimas_atomos, serial_id, destino_plots)
                
            print(f" [{status_str}]: {serial_id} | Mínima: {min_dist_global:.3f}Å | Choques: {num_choques} | Coinciden Métodos?: {coinciden_vecinos_global}")

    # --- GENERACIÓN DEL PLOT GLOBAL ---
    generar_histograma_global(distancias_minimas_globales, PLOTS_DIR)

    print(f"\n--- Experimento {FECHA_EXPERIMENTO} Finalizado con Éxito ---")
    print(f"[➔] Datos estructurados en: {EXP_OUTPUT_DIR}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Procesador Molecular con Auditoría de Distancias Mínimas.")
    parser.add_argument('rutas', nargs='*', help='Rutas de carpetas con archivos .sumviz')
    args = parser.parse_args()
    
    permitir_opcionales = preguntar_generar_reportes_y_graficos()
    
    if not args.rutas:
        ruta_por_defecto = os.path.join(BASE_DIR, "data", "input")
        print(f"[*] Usando ruta por defecto: {ruta_por_defecto}")
        ejecutar_procesamiento([ruta_por_defecto], permitir_opcionales)
    else:
        ejecutar_procesamiento(args.rutas, permitir_opcionales)