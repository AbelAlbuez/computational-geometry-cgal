#!/bin/bash
# Compila latex/informe.tex -> informe.pdf (3 pasadas para refs/citas).
# Requiere LaTeX. Si no hay pdflatex en PATH se intenta TinyTeX en $HOME
# (instalable sin sudo: wget -qO- https://yihui.org/tinytex/install-bin-unix.sh | sh).
set -e
LATEX_DIR="$(cd "$(dirname "$0")/../latex" && pwd)"

if ! command -v pdflatex >/dev/null 2>&1; then
  TT="$(echo "$HOME"/.TinyTeX/bin/* 2>/dev/null)"
  [ -d "$TT" ] && export PATH="$TT:$PATH"
fi
if ! command -v pdflatex >/dev/null 2>&1; then
  echo "No se encontró pdflatex. Instale TeX Live o TinyTeX." >&2
  exit 1
fi

cd "$LATEX_DIR"
for i in 1 2 3; do
  pdflatex -interaction=nonstopmode informe.tex >/dev/null || true
done
pdflatex -interaction=nonstopmode informe.tex | tail -3
ls -la informe.pdf
