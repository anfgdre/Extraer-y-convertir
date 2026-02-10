import os
from openbabel import openbabel

def extraer_y_convertir_quimica(nombre_archivo):
    coordenadas = []
    en_seccion = False
    
    if not os.path.exists(nombre_archivo):
        print(f"ERROR: No se encuentra el archivo en {nombre_archivo}")
        return None, None

    with open(nombre_archivo, 'r') as f:
        for linea in f:
            linea_limpia = linea.strip()
            
            # 1. Buscamos la sección 
            if "Nuclear Charges and Cartesian Coordinates" in linea:
                en_seccion = True
                continue
            
            if en_seccion:
                # Si la línea tiene datos (generalmente empieza con el símbolo del átomo)
                partes = linea_limpia.split()
                
                # Un dato válido suele tener: Atomo, Carga, X, Y, Z (5 columnas)
                if len(partes) >= 5:
                    try:
                        # Extraer símbolo químico 
                        simbolo = "".join([c for c in partes[0] if c.isalpha()])
                        
                        # Conversión Bohr a Angstrom
                        factor = 0.529177
                        x = float(partes[2]) * factor
                        y = float(partes[3]) * factor
                        z = float(partes[4]) * factor
                        
                        coordenadas.append(f"{simbolo} {x:10.6f} {y:10.6f} {z:10.6f}")
                    except ValueError:
                        # Si no se pueden convertir a float, saltamos la línea (encabezados)
                        continue
                
                # Si encontramos una línea vacía después de haber empezado a leer, terminamos
                elif not linea_limpia and len(coordenadas) > 0:
                    break

    if not coordenadas:
        print("No se encontraron coordenadas. Revisar el contenido del archivo .sumviz")
        return None, None

    # Crear el bloque de texto XYZ
    bloque_xyz = f"{len(coordenadas)}\nGenerado desde {os.path.basename(nombre_archivo)}\n" + "\n".join(coordenadas)
    
    # --- Configurar Open Babel ---
    obConv = openbabel.OBConversion()
    obConv.SetInAndOutFormats("xyz", "smi")
    mol = openbabel.OBMol()
    
    smiles_final = None
    if obConv.ReadString(mol, bloque_xyz):
        mol.ConnectTheDots()
        mol.PerceiveBondOrders()
        smiles_final = obConv.WriteString(mol).strip()
    
    return smiles_final, bloque_xyz

# --- CONFIGURACIÓN DE RUTAS ---
#--carpeta = ruta donde se buscara el archivo.sumviz
carpeta = "/home/andreamona/Documentos/Servicio social/OpenBabel/"
archivo_in = os.path.join(carpeta, "aimel_000001.sumviz") #nombre del archivo exacto
archivo_smi = os.path.join(carpeta, "resultado_final1.smi") #nombre del archivo que se creara (SMILES)
archivo_xyz = os.path.join(carpeta, "visualizacion1.xyz") #nombre del archivo que se crea (xyz)

smiles, xyz_data = extraer_y_convertir_quimica(archivo_in)

if smiles and xyz_data:
    with open(archivo_smi, "w") as f: f.write(smiles + "\n")
    with open(archivo_xyz, "w") as f: f.write(xyz_data)
    print(f"SMILES: {smiles}")
    print(f"Éxito: Se procesaron {len(xyz_data.splitlines()) - 2} átomos.")