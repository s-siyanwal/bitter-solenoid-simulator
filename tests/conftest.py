import os, sys, warnings
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "package"))
os.environ.setdefault("MPLBACKEND", "Agg")
warnings.filterwarnings("ignore", category=DeprecationWarning)
