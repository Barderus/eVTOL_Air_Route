import sys
from pathlib import Path


NEW_FOLDER = Path(__file__).resolve().parent.parent / "opensky"
if str(NEW_FOLDER) not in sys.path:
    sys.path.insert(0, str(NEW_FOLDER))

from export_data import main


if __name__ == "__main__":
    main()
