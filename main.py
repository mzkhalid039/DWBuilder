import os
import subprocess
import sys

def run_script(script_name):
    script_path = os.path.join(os.path.dirname(__file__), 'scripts', script_name)
    if os.path.exists(script_path):
        try:
            subprocess.run([sys.executable, script_path], check=True)
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running the script: {e}")
    else:
        print(f"Script {script_name} does not exist.")

def display_dwbuilder_design():
    design = """
    ***************************************************
    *                                                 *
    *                Welcome to DWBuilder             *
    *                                                 *
    ***************************************************
    *        Hey, you must know what you are doing.   *
    *  Otherwise, you might get wrong results.        *
    ***************************************************
    *     Core Developer: M.Z.Khalid                  *
    *     Main Contributors: S.M.Selbach              *
    *     Email: zeeshan.khalid039@gmail.com          *
    ***************************************************
    """
    print(design)

def display_script_options(scripts):
    print("Select a script to run:")
    script_keys = list(scripts.keys())
    half = (len(script_keys) + 1) // 2  # Calculate half, ensuring the second column gets the extra if odd

    for i in range(half):
        left = f"{script_keys[i]}: {scripts[script_keys[i]]}"
        right = f"{script_keys[i+half]}: {scripts[script_keys[i+half]]}" if i + half < len(script_keys) else ""
        print(f"{left:<30} {right}")

def main():
    display_dwbuilder_design()
    
    scripts = {
        '1': 'dwbuilder.py',
        '2': 'dbuilder.py',
        '3': 'hibuilder.py',
        '4': 'slab.py',
        '5': 'polarization.py',
        '6': 'supercell.py',
        '7': 'vasp2cif.py',
        '8': 'bondanalysis.py'

    }

    display_script_options(scripts)

    choice = input("Enter the number of the script to run: ").strip()

    if choice in scripts:
        run_script(scripts[choice])
    else:
        print("Invalid choice. Please select a valid script number.")

if __name__ == "__main__":
    main()
