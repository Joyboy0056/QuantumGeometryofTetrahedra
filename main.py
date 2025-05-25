#!/usr/bin/env python3
# This is a Claude-generated script, amazing!

"""
Main script to execute the quantum geometry tetrahedra analysis
"""

import sys
import os
from pathlib import Path

def main():
    """Main function to execute script_§5.py"""
    
    # Get the current script directory
    current_dir = Path(__file__).parent.absolute()
    
    # Path to the target script
    script_path = current_dir / "utils" / "script_§5.py"
    
    # Check if the script exists
    if not script_path.exists():
        print(f"Error: Script not found at {script_path}")
        print(f"Current directory: {current_dir}")
        print("Available files:")
        if (current_dir / "utils").exists():
            for file in (current_dir / "utils").iterdir():
                print(f"  - {file.name}")
        else:
            print("  utils/ directory not found")
        sys.exit(1)
    
    print(f"Executing: {script_path}")
    print("=" * 50)
    
    # Method 1: Using exec() to run the script in the current namespace
    try:
        # Add the script directory to Python path
        script_dir = str(script_path.parent)
        if script_dir not in sys.path:
            sys.path.insert(0, script_dir)
        
        # Change working directory to script location (in case of relative imports)
        old_cwd = os.getcwd()
        os.chdir(script_path.parent)
        
        # Execute the script
        with open(script_path, 'r', encoding='utf-8') as f:
            script_content = f.read()
        
        # Create a new namespace for the script
        script_globals = {
            '__file__': str(script_path),
            '__name__': '__main__'
        }
        
        exec(script_content, script_globals)
        
    except Exception as e:
        print(f"Error executing script: {e}")
        import traceback
        traceback.print_exc()
    
    finally:
        # Restore original working directory
        os.chdir(old_cwd)
    
    print("=" * 50)
    print("Script execution completed.")


def main_subprocess():
    """Alternative: Execute using subprocess for complete isolation"""
    import subprocess
    
    current_dir = Path(__file__).parent.absolute()
    script_path = current_dir / "utils" / "script_§5.py"
    
    if not script_path.exists():
        print(f"Error: Script not found at {script_path}")
        sys.exit(1)
    
    print(f"Executing: {script_path}")
    print("=" * 50)
    
    try:
        # Run the script as a subprocess
        result = subprocess.run([
            sys.executable, str(script_path)
        ], cwd=current_dir, capture_output=False, text=True)
        
        if result.returncode != 0:
            print(f"Script exited with code: {result.returncode}")
        
    except Exception as e:
        print(f"Error executing script: {e}")
    
    print("=" * 50)
    print("Script execution completed.")


def main_import():
    """Alternative: Import and run if the script has a main() function"""
    
    current_dir = Path(__file__).parent.absolute()
    
    # Add utils directory to Python path
    utils_dir = current_dir / "utils"
    if str(utils_dir) not in sys.path:
        sys.path.insert(0, str(utils_dir))
    
    try:
        # Try to import the script (assuming it's named script_section5.py or similar)
        # Note: Python modules can't have § in the name, so you might need to rename
        print("Attempting to import script...")
        
        # If your script is named with special characters, you might need:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "script_section5", 
            utils_dir / "script_§5.py"
        )
        script_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(script_module)
        
        # If the script has a main function, call it
        if hasattr(script_module, 'main'):
            print("Found main() function, executing...")
            script_module.main()
        else:
            print("No main() function found in script")
            
    except Exception as e:
        print(f"Error importing/executing script: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    print("Quantum Geometry Tetrahedra - Section 5 Analysis")
    print("=" * 50)
    
    # Choose execution method
    method = "exec"  # Options: "exec", "subprocess", "import"
    
    if method == "exec":
        main()
    elif method == "subprocess":
        main_subprocess()
    elif method == "import":
        main_import()
    else:
        print(f"Unknown execution method: {method}")
        sys.exit(1)