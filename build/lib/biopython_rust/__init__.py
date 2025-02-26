try:
    from biopython_rust import __version__
except ImportError:
    __version__ = "0.1.0"

# Import submodules
from . import seq
from . import io