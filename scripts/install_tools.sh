#!/usr/bin/env bash
set -euo pipefail

echo "[install] apt packages…"
sudo apt-get update
sudo apt-get install -y \
  fastqc hisat2 samtools subread cutadapt wget gzip gawk

echo "[install] Python tools…"
python3 -m pip install --user --upgrade pip
python3 -m pip install --user "multiqc==1.31"

# Ensure MultiQC is on PATH for future shells
if ! grep -q '.local/bin' "$HOME/.bashrc"; then
  echo 'export PATH="$HOME/.local/bin:$PATH"' >> "$HOME/.bashrc"
fi

echo "[done] Open a new terminal (or run: source ~/.bashrc)."
