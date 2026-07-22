#!/usr/bin/env python3

"""Install LMPTorch as a selectable package in a recent LAMMPS source tree."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


PACKAGE_NAME = "LMPTORCH"
SOURCE_NAMES = ("fix_python_torch.cpp", "fix_python_torch.h")


def insert_once(text: str, marker: str, addition: str, description: str) -> str:
    if addition in text:
        return text
    if marker not in text:
        raise RuntimeError(f"Could not locate {description} marker: {marker!r}")
    return text.replace(marker, marker + addition, 1)


def install_package(lammps_dir: Path) -> None:
    repo_dir = Path(__file__).resolve().parent
    package_source = repo_dir / "src" / PACKAGE_NAME
    module_source = repo_dir / "cmake" / "Modules" / "Packages" / f"{PACKAGE_NAME}.cmake"
    cmake_lists = lammps_dir / "cmake" / "CMakeLists.txt"
    lammps_src = lammps_dir / "src"

    if not cmake_lists.is_file() or not (lammps_src / "version.h").is_file():
        raise RuntimeError(f"Not a LAMMPS source tree: {lammps_dir}")

    package_target = lammps_src / PACKAGE_NAME
    package_target.mkdir(parents=True, exist_ok=True)
    for name in SOURCE_NAMES:
        source = package_source / name
        target = package_target / name
        shutil.copy2(source, target)

        # Remove files left by the legacy copy-to-src installation method.
        legacy_target = lammps_src / name
        if legacy_target.exists():
            if legacy_target.read_bytes() != source.read_bytes():
                raise RuntimeError(f"Refusing to replace modified legacy file: {legacy_target}")
            legacy_target.unlink()

    module_target = lammps_dir / "cmake" / "Modules" / "Packages" / module_source.name
    shutil.copy2(module_source, module_target)

    text = cmake_lists.read_text(encoding="utf-8")
    text = insert_once(text, "  LEPTON\n", "  LMPTORCH\n", "standard package list")
    if "PYTHON LMPTORCH ML-IAP" not in text:
        marker = "PYTHON ML-IAP"
        if marker not in text:
            raise RuntimeError("Could not locate the package-module include list")
        text = text.replace(marker, "PYTHON LMPTORCH ML-IAP", 1)
    cmake_lists.write_text(text, encoding="utf-8")

    print(f"Installed {PACKAGE_NAME} package into {lammps_dir}")
    print("Configure with PKG_PYTHON=ON and PKG_LMPTORCH=ON")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("lammps_dir", type=Path, help="LAMMPS source tree")
    args = parser.parse_args()
    install_package(args.lammps_dir.resolve())


if __name__ == "__main__":
    main()
