#!/bin/bash
# Descobre o diretório onde está o script
FORMATTER_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Caminho para o Julia (ajuste se não estiver no PATH)
JULIA_BIN="julia"

# Executa o format.jl com o projeto local
$JULIA_BIN --project="$FORMATTER_DIR" "$FORMATTER_DIR/format.jl" "$@"
