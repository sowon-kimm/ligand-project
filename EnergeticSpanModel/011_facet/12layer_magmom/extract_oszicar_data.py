import os
import glob
import pandas as pd

def extract_oszicar_data_with_pattern():
    results = []   
    search_pattern = os.path.join('.', '*', 'OSZICAR')
    file_list = glob.glob(search_pattern)

    if not file_list:
        print("path error")
        return

    for file_path in file_list:
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                
                if not lines:
                    continue

                last_line = lines[-1]
                parts = last_line.split()
                
                if "E0=" in parts:
                    idx = parts.index("E0=") + 1
                    e0_val = parts[idx]                   
                    folder_name = os.path.dirname(file_path)
                    
                    results.append({
                        "path": folder_name,
                        "e0": float(e0_val)
                    })
                    print(f"[success] {folder_name} : {e0_val}")
                    
        except Exception as e:
            print(f"Error reading {file_path}: {e}")

    if results:
        df = pd.DataFrame(results)
        df = df.sort_values(by="e0")
        output_file = "VASP_OSZICAR_Top_Summary.xlsx"
        df.to_excel(output_file, index=False)
        
        print("\n--- summary ---")
        print(df.head())
        print(f"\nsaved in {output_file}.")
    else:
        print("error")

if __name__ == "__main__":
    extract_oszicar_data_with_pattern()
