import os
from openbabel import openbabel

# Configuración dinámica de rutas
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUT_DIR = os.path.join(BASE_DIR, "data", "input")
OUTPUT_DIR = os.path.join(BASE_DIR, "data", "output")

def extraer_y_convertir_quimica(nombre_archivo):
    coordenadas = []
    en_seccion = False
    factor_bohr_to_angstrom = 0.529177
    
    if not os.path.exists(nombre_archivo):
        raise FileNotFoundError(f"No se encontró: {nombre_archivo}")

    with open(nombre_archivo, 'r') as f:
        for linea in f:
            linea_limpia = linea.strip()
            
            if "Nuclear Charges and Cartesian Coordinates" in linea:
                en_seccion = True
                continue
            
            if en_seccion:
                partes = linea_limpia.split()
                if len(partes) >= 5:
                    try:
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

    bloque_xyz = f"{len(coordenadas)}\nGenerado desde {os.path.basename(nombre_archivo)}\n" + "\n".join(coordenadas)
    
    # --- Open Babel Logic ---
    obConv = openbabel.OBConversion()
    obConv.SetInAndOutFormats("xyz", "smi")
    mol = openbabel.OBMol()
    
    smiles_final = None
    if obConv.ReadString(mol, bloque_xyz):
        mol.ConnectTheDots()
        mol.PerceiveBondOrders()
        smiles_final = obConv.WriteString(mol).strip()
    
    return smiles_final, bloque_xyz

def ejecutar_procesamiento(archivo_nombre):
    # Asegurar que la carpeta de salida exista
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    ruta_in = os.path.join(INPUT_DIR, archivo_nombre)
    nombre_base = os.path.splitext(archivo_nombre)[0]
    
    try:
        smiles, xyz_data = extraer_y_convertir_quimica(ruta_in)
        if smiles:
            with open(os.path.join(OUTPUT_DIR, f"{nombre_base}.smi"), "w") as f: f.write(smiles + "\n")
            with open(os.path.join(OUTPUT_DIR, f"{nombre_base}.xyz"), "w") as f: f.write(xyz_data)
            print(f"✅ Procesado con éxito: {archivo_nombre}")
    except Exception as e:
        print(f"❌ Error en {archivo_nombre}: {e}")

if __name__ == "__main__":
    # Ejemplo: Procesar un archivo específico
    ejecutar_procesamiento("aimel_000001.sumviz")