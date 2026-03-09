import os
import glob
from openbabel import openbabel

# --- CONFIGURACIÓN DINÁMICA DE RUTAS ---
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(_file_)))
INPUT_DIR = os.path.join(BASE_DIR, "data", "input")
OUTPUT_DIR = os.path.join(BASE_DIR, "data", "output")

def extraer_coordenadas_sumviz(ruta_archivo):
    """
    Función independiente: Extrae puramente las coordenadas de un .sumviz.
    Retorna un string formateado en XYZ o None si falla.
    """
    coordenadas = []
    en_seccion = False
    factor_bohr_to_angstrom = 0.529177
    
    try:
        with open(ruta_archivo, 'r') as f:
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
            return None

        # Retorna el bloque XYZ listo para escribir o procesar
        return f"{len(coordenadas)}\nExtraído de {os.path.basename(ruta_archivo)}\n" + "\n".join(coordenadas)
    except Exception as e:
        print(f"Error en extracción: {e}")
        return None

def convertir_xyz_a_smiles(bloque_xyz):
    """
    Función independiente: Recibe un string XYZ y devuelve un SMILES usando Open Babel.
    """
    try:
        obConv = openbabel.OBConversion()
        obConv.SetInAndOutFormats("xyz", "smi")
        mol = openbabel.OBMol()
        
        if obConv.ReadString(mol, bloque_xyz):
            mol.ConnectTheDots()
            mol.PerceiveBondOrders()
            return obConv.WriteString(mol).strip()
    except Exception as e:
        print(f"Error en Open Babel: {e}")
    return None

def ejecutar_procesamiento(incluir_smiles=False):
    """
    Orquestador: Recorre la carpeta y decide qué funciones llamar.
    """
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    archivos = glob.glob(os.path.join(INPUT_DIR, "*.sumviz"))
    
    if not archivos:
        print(f"\n[!] No hay archivos .sumviz en {INPUT_DIR}")
        return

    print(f"\n--- Procesando {len(archivos)} archivos ---")
    for ruta in archivos:
        nombre = os.path.splitext(os.path.basename(ruta))[0]
        
        # 1. Siempre extraemos XYZ
        xyz_data = extraer_coordenadas_sumviz(ruta)
        
        if xyz_data:
            # Guardar XYZ
            with open(os.path.join(OUTPUT_DIR, f"{nombre}.xyz"), "w") as f:
                f.write(xyz_data)
            
            # 2. Solo si el usuario pidió SMILES
            if incluir_smiles:
                smiles = convertir_xyz_a_smiles(xyz_data)
                if smiles:
                    with open(os.path.join(OUTPUT_DIR, f"{nombre}.smi"), "w") as f:
                        f.write(smiles + "\n")
                    print(f"[OK] {nombre}: XYZ + SMILES")
                else:
                    print(f"[!] {nombre}: XYZ guardado, pero falló SMILES")
            else:
                print(f"[OK] {nombre}: Solo XYZ")
        else:
            print(f"[ERROR] {nombre}: No se detectaron coordenadas")

def menu():
    while True:
        print("\n" + "="*30)
        print("  EXTRACTOR DE COORDENADAS  ")
        print("="*30)
        print("1. Extraer solo archivos XYZ")
        print("2. Extraer XYZ y generar SMILES")
        print("3. Salir")
        opcion = input("\nSelecciona una opción: ")

        if opcion == "1":
            ejecutar_procesamiento(incluir_smiles=False)
        elif opcion == "2":
            ejecutar_procesamiento(incluir_smiles=True)
        elif opcion == "3":
            print("Saliendo...")
            break
        else:
            print("Opción no válida.")

if _name_ == "_main_":
    menu()