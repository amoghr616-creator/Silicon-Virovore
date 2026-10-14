"""Independent target-contact, coordinate, and clash evaluation for a pose."""

try:
    from .target_audit import main
except ImportError:  # direct script execution
    from target_audit import main


if __name__ == "__main__":
    raise SystemExit(main())
