import os
import json
import pandas as pd
import glob
from collections import Counter #

# --- CONFIGURACIÓN ---
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_BASE_DIR = os.path.join(BASE_DIR, "data", "output")
ARCHIVO_MAESTRO_CSV = os.path.join(OUTPUT_BASE_DIR, "Reporte_Maestro_Moleculas.csv")

def consolidar_datos():
    todos_los_datos = []
    total_anomalias_detectadas = 0
    censo_global_atomos = Counter() # Para el conteo total de tipos y cantidades
    
    patron_busqueda = os.path.join(OUTPUT_BASE_DIR, "Ejecucion_*", "auditoria_final.json")
    archivos_json = glob.glob(patron_busqueda)
    
    if not archivos_json:
        print("No se encontraron archivos de auditoría para consolidar.")
        return

    print(f"Encontrados {len(archivos_json)} reportes. Procesando...")

    for ruta_json in archivos_json:
        nombre_sesion = os.path.basename(os.path.dirname(ruta_json))
        try:
            with open(ruta_json, 'r', encoding='utf-8') as f:
                datos_sesion = json.load(f)
                for molecula in datos_sesion:
                    molecula['sesion_id'] = nombre_sesion
                    total_anomalias_detectadas += molecula.get('anomalias_vecinos', 0)
                    
                    # --- CONTEO GLOBAL DE ÁTOMOS ---
                    formula_dict = molecula.get('formula', {})
                    censo_global_atomos.update(formula_dict) # Suma tipos y cantidades al total
                    
                    molecula['formula_txt'] = "".join([f"{k}{v}" for k, v in formula_dict.items()])
                    todos_los_datos.append(molecula)
        except Exception as e:
            print(f"Error procesando {ruta_json}: {e}")

    if todos_los_datos:
        df = pd.DataFrame(todos_los_datos)
        columnas = ['sesion_id', 'archivo', 'status', 'formula_txt', 'distancia_max', 'anomalias_vecinos', 'smiles_ob', 'smiles_rd']
        df = df[columnas]
        df.to_csv(ARCHIVO_MAESTRO_CSV, index=False, encoding='utf-8-sig')
        
        print(f"\nReporte consolidado: {ARCHIVO_MAESTRO_CSV}")
        print("-" * 40)
        print(f"Total de moléculas: {len(df)}")
        print(f"Total de anomalías: {total_anomalias_detectadas}")
        print("\nCENSO GLOBAL DE ÁTOMOS PROCESADOS:") #
        for elemento, cantidad in sorted(censo_global_atomos.items()):
            print(f"  - {elemento}: {cantidad} átomos") # Muestra qué tipos y cuántos son
        print("-" * 40)

if __name__ == "__main__":
    consolidar_datos()