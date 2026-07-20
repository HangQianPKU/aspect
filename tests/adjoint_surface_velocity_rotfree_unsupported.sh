#!/bin/bash

set -e

output_file="$(mktemp)"
cat > "${output_file}"

grep -q "net-rotation projection is not implemented" "${output_file}"
echo "unsupported rot-free surface velocity objective failed with the expected message"
