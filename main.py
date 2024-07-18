# main.py
import os
import subprocess
import sys

def run_script(script_name):
    if os.path.exists(script_name):
        try:
            subprocess.run([sys.executable, script_name], check=True)
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running the script: {e}")
    else:
        print(f"Script {script_name} does not exist.")

def main():
    scripts = {
        '1': 'dwbuilder.py',
        '2': 'dbuilder.py',
        '3': 'hibuilder.py',
        '4': 'slab.py'
    }

    print("Select a script to run:")
    for key, value in scripts.items():
        print(f"{key}: {value}")

    choice = input("Enter the number of the script to run: ").strip()

    if choice in scripts:
        run_script(scripts[choice])
    else:
        print("Invalid choice. Please select a valid script number.")

if __name__ == "__main__":
    main()
