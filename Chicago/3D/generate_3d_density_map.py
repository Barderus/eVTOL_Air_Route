import sys
from pathlib import Path


NEW_FOLDER = Path(__file__).resolve().parent
if str(NEW_FOLDER) not in sys.path:
    sys.path.insert(0, str(NEW_FOLDER))

# Keep the historical entry point while delegating the real work to the main
# generator script next to this wrapper.
from make_3d_density_map import main


if __name__ == "__main__":
    main()
