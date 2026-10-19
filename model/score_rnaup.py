"""Backward-compatible entry point → ``model.energy.score_rnaup``."""

from model.energy.score_rnaup import *  # noqa: F403
from model.energy.score_rnaup import main

if __name__ == "__main__":
    main()
