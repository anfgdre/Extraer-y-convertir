import os
import glob
from openbabel import openbabel

# --- CONFIGURACIÓN DINÁMICA DE RUTAS ---
# Localiza la raíz del proyecto respecto a este archivo
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUT_DIR = os.path.join(BASE_DIR, "data", "input")
OUTPUT_DIR = os.path.join(BASE_DIR, "data", "output")

def extraer_y_convertir_quimica(ruta_archivo):
    coordenadas = []
    en_seccion = False
    factor_bohr_to_angstrom = 0.529177
    
    try:
        with open(ruta_archivo, 'r') as f:
            for linea in f:
                linea_limpia = linea.strip()
                
                # Identificar sección de coordenadas
                if "Nuclear Charges and Cartesian Coordinates" in linea:
                    en_seccion = True
                    continue
                
                if en_seccion:
                    partes = linea_limpia.split()
                    if len(partes) >= 5:
                        try:
                            # Extraer símbolo (letras) y convertir unidades
                            simbolo = "".join([c for c in partes[0] if c.isalpha()])
                            x = float(partes[2]) * factor_bohr_to_angstrom
                            y = float(partes[3]) * factor_bohr_to_angstrom
                            z = float(partes[4]) * factor_bohr_to_angstrom
                            coordenadas.append(f"{simbolo} {x:10.6f} {y:10.6f} {z:10.6f}")
                        except ValueError:
                            continue
                    elif not linea_limpia and len(coordenadas) > 0:
                        break

        if not coordenadas:
            return None, None

        # Formato XYZ
        bloque_xyz = f"{len(coordenadas)}\nGenerado desde {os.path.basename(ruta_archivo)}\n" + "\n".join(coordenadas)
        
        # --- Lógica de Open Babel ---
        obConv = openbabel.OBConversion()
        obConv.SetInAndOutFormats("xyz", "smi")
        mol = openbabel.OBMol()
        
        smiles_final = None
        if obConv.ReadString(mol, bloque_xyz):
            mol.ConnectTheDots()
            mol.PerceiveBondOrders()
            smiles_final = obConv.WriteString(mol).strip()
            
        return smiles_final, bloque_xyz

    except Exception as e:
        print(f"Error procesando el archivo: {e}")
        return None, None

def procesar_todo():
    # Asegurar que existan las carpetas
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Buscar todos los archivos .sumviz en la carpeta input
    archivos_sumviz = glob.glob(os.path.join(INPUT_DIR, "*.sumviz"))
    
    if not archivos_sumviz:
        print(f"No se encontraron archivos .sumviz en: {INPUT_DIR}")
        return

    print(f"--- Iniciando procesamiento de {len(archivos_sumviz)} archivos ---")
    
    for ruta_in in archivos_sumviz:
        nombre_base = os.path.splitext(os.path.basename(ruta_in))[0]
        smiles, xyz_data = extraer_y_convertir_quimica(ruta_in)
        
        if smiles and xyz_data:
            # Guardar resultados
            with open(os.path.join(OUTPUT_DIR, f"{nombre_base}.smi"), "w") as f:
                f.write(smiles + "\n")
            with open(os.path.join(OUTPUT_DIR, f"{nombre_base}.xyz"), "w") as f:
                f.write(xyz_data)
            print(f"Éxito: {nombre_base}")
        else:
            print(f"Falló: {nombre_base} (No se extrajeron coordenadas)")

if __name__ == "__main__":
    procesar_todo()